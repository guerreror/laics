/*
 *  smc.cpp
 *  SMC entry point (minimal stub).
 *
 *  Runs the standard setup and ARG simulation, then emits the first
 *  site tree and exits. 
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <functional>
#include <chrono>
#include <random>
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <cstdint>
#include <cstring>
#include <unordered_map>

using namespace std;

#include "parameters.h"
#include "world.h"
#include "ran_mk.h"
#include "sitenode.h"
#include "snptree.h"
#include "migprob.h"
#include "treemod.h"
#include "simulate_smc.h"
#include "chromosome.h"
#include "smc_helpers.h"

std::random_device rd;
auto seed = rd();
std::mt19937_64 gen(seed);

static bool looksLikePathArg(const string& s)
{
    return s == "." || s.find('/') != string::npos || s.find('\\') != string::npos;
}

static string pathJoin(const string& dir, const string& file)
{
    if (dir.empty() || dir == ".") return file;
    if (dir.back() == '/' || dir.back() == '\\') return dir + file;
    return dir + "/" + file;
}

static vector<vector<double>> buildMigMatrix(Parameters &p)
{
    unsigned int nPops = static_cast<int>(p.paramData->popSizeVec.size());
    vector<vector<double>> mig_prob;
    mig_prob.resize(nPops);
    double m = p.paramData->migRate.at(0) / (2 * p.paramData->totalPopSize);

    for (unsigned int i = 0; i < nPops; i++)
        mig_prob.at(i).resize(nPops);
    for (unsigned int i = 0; i < nPops; i++)
    {
        for (unsigned int j = 0; j < nPops; j++)
        {
            if (nPops == 1)
                mig_prob.at(i).at(j) = 1;
            else if (i == j)
                mig_prob.at(i).at(j) = 1 - m;
            else if (i == j - 1)
            {
                mig_prob.at(i).at(j) = m / 2;
                if (i == 0)
                    mig_prob.at(i).at(j) = m;
            }
            else if (i == j + 1)
            {
                mig_prob.at(i).at(j) = m / 2;
                if (i == nPops - 1)
                    mig_prob.at(i).at(j) = m;
            }
            else
                mig_prob.at(i).at(j) = 0;
        }
    }
    return mig_prob;
}

static vector<vector<double>> adaptMatrixForPops(const vector<vector<double>> &m, unsigned int nPops)
{
    if (m.size() == nPops) {
        return m;
    }
    // Some migration schedules include an ancestral/root deme at index 0.
    // For simulations with only extant populations, strip row/col 0.
    if (m.size() == nPops + 1) {
        bool square = true;
        for (const auto &row : m) {
            if (row.size() != m.size()) {
                square = false;
                break;
            }
        }
        if (square) {
            vector<vector<double>> out(nPops, vector<double>(nPops, 0.0));
            for (unsigned int i = 0; i < nPops; ++i) {
                for (unsigned int j = 0; j < nPops; ++j) {
                    out[i][j] = m[i + 1][j + 1];
                }
            }
            return out;
        }
    }
    return m;
}
// Check if any migration events are present in the schedule
static bool scheduleHasMigration(
    const vector<pair<double, Matrix>>& schedule)
{
    for (const auto& item : schedule) {
        const Matrix& matrix = item.second;
        for (size_t i = 0; i < matrix.size(); ++i) {
            for (size_t j = 0; j < matrix[i].size(); ++j) {
                // If any off-diagonal entry is greater than zero, migration is present
                if (i != j && matrix[i][j] > 0.0) {
                    return true;
                }
            }
        }
    }
    return false;
}

static double sampleHopDeltaFromRho(double rho)
{
    if (rho <= 0.0) return -1.0;
    return randexp(rho);
}

static bool pickWeightedEdge(const vector<EdgeWeight>& standard_edges,
                             const vector<EdgeWeight>& inverted_edges,
                             unsigned long* outParent,
                             unsigned long* outChild)
{
    if (!outParent || !outChild) return false;
    vector<EdgeWeight> all;
    all.reserve(standard_edges.size() + inverted_edges.size());
    for (const auto& e : standard_edges) if (e.weight > 0.0) all.push_back(e);
    for (const auto& e : inverted_edges) if (e.weight > 0.0) all.push_back(e);
    if (all.empty()) return false;

    double totalW = 0.0;
    for (const auto& e : all) totalW += e.weight;
    if (totalW <= 0.0) return false;

    const double roll = randreal(0, totalW);
    double cum = 0.0;
    for (const auto& e : all) {
        cum += e.weight;
        if (roll <= cum) {
            *outParent = e.parent;
            *outChild = e.child;
            return true;
        }
    }
    *outParent = all.back().parent;
    *outChild = all.back().child;
    return true;
}

static void collectEdgeWeightsFromTree(
    TreeNode* node,
    const Parameters::ParameterData& params,
    double r,
    vector<EdgeWeight>& standard_edges,
    vector<EdgeWeight>& inverted_edges)
{
    if (!node) return;
    for (auto* child : node->children) {
        if (!child) continue;

        // Cutting an edge on a unary root stem removes the entire genealogy,
        // leaving no retained tree for the lineage to reattach to.
        TreeNode* branchingAncestor = node;
        while (branchingAncestor && branchingAncestor->children.size() == 1) {
            branchingAncestor = branchingAncestor->parent;
        }
        const bool leavesRetainedTree =
            branchingAncestor && branchingAncestor->children.size() > 1;

        const unsigned int pop = child->context.pop;
        const SMCActiveState state = activeStateAtTime_SMC(params, child->time);
        if (leavesRetainedTree && pop < state.popSizes.size() && pop < state.invFreqs.size()) {
            const double pI = state.invFreqs[pop];
            const double branchL = node->time - child->time;
            if (branchL > 0.0) {
                const double base = r * branchL;
                if (child->context.inversion == 1) {
                    inverted_edges.push_back({node->id, child->id, base * pI});
                } else {
                    standard_edges.push_back({node->id, child->id, base * (1.0 - pI)});
                }
            }
        }
        collectEdgeWeightsFromTree(child, params, r, standard_edges, inverted_edges);
    }
}

static void writeEdgeWeightsCSV(const vector<EdgeWeight>& edges, const std::string& path)
{
    std::ofstream out(path.c_str());
    if (!out.is_open()) return;
    out << "parent,child,weight\n";
    for (const auto& e : edges) {
        out << e.parent << "," << e.child << "," << e.weight << "\n";
    }
}

static void writeTreeArtifacts(TreeNode* tree, const std::string& base)
{
    const std::string dot = base + ".dot";
    const std::string collapsedDot = base + "_collapsed.dot";
    writeTreeDOT(tree, dot);
    writeCollapsedTreeDOT(tree, collapsedDot);

    std::ostringstream cmd;
    cmd << "python3 tools/dot_to_tskit_png.py "
        << dot
        << " --svg " << base << ".tskit.svg";
    const int rc = std::system(cmd.str().c_str());
    if (rc != 0) {
        std::cerr << "Warning: failed to render " << dot
                  << " to tskit SVG.\n";
    }

    std::ostringstream cmdCollapsed;
    cmdCollapsed << "python3 tools/dot_to_tskit_png.py "
                 << collapsedDot
                 << " --svg " << base << "_collapsed.tskit.svg";
    const int rcCollapsed = std::system(cmdCollapsed.str().c_str());
    if (rcCollapsed != 0) {
        std::cerr << "Warning: failed to render " << collapsedDot
                  << " to tskit SVG.\n";
    }
}

static void writeTreeSnapshotCSVRows(std::ofstream& out,
                                     TreeNode* node,
                                     int run,
                                     int hop,
                                     double xStart,
                                     double xEnd,
                                     long long parentId)
{
    if (!node) return;
    out << run << ","
        << hop << ","
        << std::setprecision(17) << xStart << ","
        << std::setprecision(17) << xEnd << ","
        << node->id << ","
        << parentId << ","
        << std::setprecision(17) << node->time << ","
        << node->context.pop << ","
        << node->context.inversion << "\n";
    for (auto* child : node->children) {
        writeTreeSnapshotCSVRows(out, child, run, hop, xStart, xEnd,
                                 static_cast<long long>(node->id));
    }
}

static void appendTreeSnapshotCSV(std::ofstream& out,
                                  TreeNode* tree,
                                  int run,
                                  int hop,
                                  double xStart,
                                  double xEnd)
{
    if (!out.is_open()) return;
    writeTreeSnapshotCSVRows(out, tree, run, hop, xStart, xEnd, -1);
}

template <typename T>
static void packBinaryValue(char* buffer, size_t& offset, const T& value)
{
    std::memcpy(buffer + offset, &value, sizeof(T));
    offset += sizeof(T);
}

static void writeTreeSnapshotBinaryRows(std::ofstream& out,
                                        TreeNode* node,
                                        int32_t run,
                                        int32_t hop,
                                        double xStart,
                                        double xEnd,
                                        int64_t parentId)
{
    if (!node) return;
    const uint64_t nodeId = static_cast<uint64_t>(node->id);
    const double time = node->time;
    const uint32_t pop = static_cast<uint32_t>(node->context.pop);
    const uint16_t inversion = static_cast<uint16_t>(node->context.inversion);
    char record[sizeof(run) + sizeof(hop) + sizeof(xStart) + sizeof(xEnd) +
                sizeof(nodeId) + sizeof(parentId) + sizeof(time) +
                sizeof(pop) + sizeof(inversion)];
    size_t offset = 0;
    packBinaryValue(record, offset, run);
    packBinaryValue(record, offset, hop);
    packBinaryValue(record, offset, xStart);
    packBinaryValue(record, offset, xEnd);
    packBinaryValue(record, offset, nodeId);
    packBinaryValue(record, offset, parentId);
    packBinaryValue(record, offset, time);
    packBinaryValue(record, offset, pop);
    packBinaryValue(record, offset, inversion);
    out.write(record, sizeof(record));
    for (auto* child : node->children) {
        writeTreeSnapshotBinaryRows(out, child, run, hop, xStart, xEnd,
                                    static_cast<int64_t>(node->id));
    }
}

static void appendTreeSnapshotBinary(std::ofstream& out,
                                     TreeNode* tree,
                                     int run,
                                     int hop,
                                     double xStart,
                                     double xEnd)
{
    if (!out.is_open()) return;
    writeTreeSnapshotBinaryRows(out, tree, static_cast<int32_t>(run),
                                static_cast<int32_t>(hop), xStart, xEnd, -1);
}

static std::string formatCoordForFilename(double x)
{
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(6) << x;
    std::string out = ss.str();
    while (!out.empty() && out.back() == '0') out.pop_back();
    if (!out.empty() && out.back() == '.') out.pop_back();
    std::replace(out.begin(), out.end(), '.', 'p');
    std::replace(out.begin(), out.end(), '-', 'm');
    if (out.empty()) out = "0";
    return out;
}

struct TreeShapeStats {
    unsigned long node_count = 0;
    unsigned long edge_count = 0;
    unsigned long leaf_count = 0;
    unsigned long unary_node_count = 0;
    unsigned long branching_node_count = 0;
    unsigned long max_depth = 0;
    double root_time = 0.0;
    double max_time = 0.0;
    double total_branch_length = 0.0;
    double unary_parent_branch_length = 0.0;
    double terminal_leaf_branch_length = 0.0;
    double internal_branching_branch_length = 0.0;
    double terminal_from_unary_parent_branch_length = 0.0;
    double terminal_from_branching_parent_branch_length = 0.0;
    double internal_unary_chain_branch_length = 0.0;
    double internal_branching_exclusive_branch_length = 0.0;
    double standard_branch_length = 0.0;
    double inverted_branch_length = 0.0;
};

static void accumulateTreeShapeStats(TreeNode* node,
                                     unsigned long depth,
                                     TreeShapeStats& stats)
{
    if (!node) return;

    stats.node_count++;
    stats.max_depth = std::max(stats.max_depth, depth);
    stats.max_time = std::max(stats.max_time, node->time);

    if (node->children.empty()) {
        stats.leaf_count++;
    } else if (node->children.size() == 1) {
        stats.unary_node_count++;
    } else {
        stats.branching_node_count++;
    }

    for (auto* child : node->children) {
        if (!child) continue;
        const double branchL = node->time - child->time;
        if (branchL > 0.0) {
            stats.edge_count++;
            stats.total_branch_length += branchL;
            if (child->context.inversion == 1) {
                stats.inverted_branch_length += branchL;
            } else {
                stats.standard_branch_length += branchL;
            }

            // Diagnostic only. The first few summaries intentionally overlap
            // historical questions we were asking; the exclusive categories
            // below classify each branch exactly once for clean composition plots.
            if (node->children.size() == 1) {
                stats.unary_parent_branch_length += branchL;
            }
            if (child->children.empty()) {
                stats.terminal_leaf_branch_length += branchL;
            } else if (child->children.size() > 1) {
                stats.internal_branching_branch_length += branchL;
            }

            if (child->children.empty() && node->children.size() == 1) {
                stats.terminal_from_unary_parent_branch_length += branchL;
            } else if (child->children.empty()) {
                stats.terminal_from_branching_parent_branch_length += branchL;
            } else if (node->children.size() == 1) {
                stats.internal_unary_chain_branch_length += branchL;
            } else {
                stats.internal_branching_exclusive_branch_length += branchL;
            }
        }
        accumulateTreeShapeStats(child, depth + 1, stats);
    }
}

static TreeShapeStats summarizeTreeShape(TreeNode* root)
{
    TreeShapeStats stats;
    if (root) {
        stats.root_time = root->time;
        accumulateTreeShapeStats(root, 0, stats);
    }
    return stats;
}

static void writeTreeShapeHeader(const string& path)
{
    std::ofstream out(path.c_str());
    if (!out.is_open()) return;
    out << "run,hop,current_x,rho,root_time,max_time,node_count,edge_count,"
        << "leaf_count,unary_node_count,branching_node_count,max_depth,"
        << "total_branch_length,unary_parent_branch_length,"
        << "terminal_leaf_branch_length,internal_branching_branch_length,"
        << "terminal_from_unary_parent_branch_length,"
        << "terminal_from_branching_parent_branch_length,"
        << "internal_unary_chain_branch_length,"
        << "internal_branching_exclusive_branch_length,"
        << "standard_branch_length,inverted_branch_length\n";
}

static void appendTreeShapeRow(const string& path,
                               int run,
                               int hop,
                               double currentX,
                               double rho,
                               const TreeShapeStats& stats)
{
    std::ofstream out(path.c_str(), std::ios::app);
    if (!out.is_open()) return;
    out << run << ","
        << hop << ","
        << currentX << ","
        << rho << ","
        << stats.root_time << ","
        << stats.max_time << ","
        << stats.node_count << ","
        << stats.edge_count << ","
        << stats.leaf_count << ","
        << stats.unary_node_count << ","
        << stats.branching_node_count << ","
        << stats.max_depth << ","
        << stats.total_branch_length << ","
        << stats.unary_parent_branch_length << ","
        << stats.terminal_leaf_branch_length << ","
        << stats.internal_branching_branch_length << ","
        << stats.terminal_from_unary_parent_branch_length << ","
        << stats.terminal_from_branching_parent_branch_length << ","
        << stats.internal_unary_chain_branch_length << ","
        << stats.internal_branching_exclusive_branch_length << ","
        << stats.standard_branch_length << ","
        << stats.inverted_branch_length << "\n";
}

static void printProgress(int completed, int total, int& lastBucket)
{
    if (total <= 0) return;
    int bucket = (completed >= total) ? 10 : static_cast<int>((static_cast<long long>(completed) * 10) / total);
    if (bucket == lastBucket) return;
    lastBucket = bucket;

    const int width = 40;
    const double frac = static_cast<double>(completed) / static_cast<double>(total);
    int filled = static_cast<int>(frac * width);
    if (filled > width) filled = width;

    std::cerr << "\rProgress: [";
    for (int i = 0; i < width; ++i) {
        std::cerr << (i < filled ? '#' : '-');
    }
    std::cerr << "] " << static_cast<int>(frac * 100.0)
              << "% (" << completed << "/" << total << " replicates)";
    if (completed >= total) {
        std::cerr << "\n";
    }
    std::cerr.flush();
}

int main(int argc, const char *argv[])
{
    std::cerr << "Random Seed: " << seed << '\n';

    auto input_seed = std::stoi(argv[1]);
    if (input_seed != -1)
    {
        std::cerr << "Input Seed: " << input_seed << '\n';
        seed = input_seed;
        gen.seed(input_seed);
    }

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    stringstream infile;
    infile << "inLABP.pars";

    int param_arg_end = argc;
    string output_dir = ".";
    if (argc > 2 && looksLikePathArg(argv[argc - 1])) {
        output_dir = argv[argc - 1];
        param_arg_end = argc - 1;
    }
    const string tree_dir = (output_dir == ".") ? "." : pathJoin(output_dir, "trees");

    vector<string> param_vec(argv + 2, argv + param_arg_end);
    Parameters params(infile.str().c_str(), param_vec);

    unsigned int nRuns = params.paramData->nRuns;
    unsigned int nSites = params.paramData->n_SNPs;
    const bool targetMode = !params.paramData->targetSNPs.empty();
    const bool writeAllDiagnostics = !targetMode && params.paramData->smcVerbose;
    const string tree_shape_path = pathJoin(output_dir, "smc_tree_shape.csv");
    const string coalescence_diagnostic_path =
        pathJoin(output_dir, "smc_coalescence_diagnostics.csv");
    bool treeShapeHeaderWritten = false;

    {
        std::ofstream coalescenceDiagnostics(coalescence_diagnostic_path);
        if (coalescenceDiagnostics.is_open()) {
            coalescenceDiagnostics
                << "run,hop,current_x,phase,lineage_time,population,arrangement,"
                << "population_size,arrangement_frequency,context_size,lineage_count,"
                << "eligible_pair_count,pair_rate_used,total_c_used,total_m,total_g\n";
        }
    }

    std::cerr << "\n\n";
    int progressBucket = -1;
    printProgress(0, static_cast<int>(nRuns), progressBucket);

    const std::string mig_json = "src/migration_matrices.json";
    const std::string adjusted_standard_mig_json =
        "src/adjusted_standard_migration_matrices.json";
    const std::string adjusted_inverted_mig_json =
        "src/adjusted_inverted_migration_matrices.json";
    auto schedule = readMigrationSchedule(mig_json);
    auto adjusted_standard_schedule =
        readMigrationSchedule(adjusted_standard_mig_json);
    auto adjusted_inverted_schedule =
        readMigrationSchedule(adjusted_inverted_mig_json);

    // Check if migration is enabled based on the schedule
    bool migration_enabled = scheduleHasMigration(schedule);
    if (schedule.empty()) {
        for (double rate : params.paramData->migRate) {
            if (rate > 0.0) {
                migration_enabled = true;
                break;
            }
        }
    }
    // Check if the required adjusted migration schedules are available
    if (migration_enabled &&
        (adjusted_standard_schedule.empty() || adjusted_inverted_schedule.empty())) {
        std::cerr
            << "Error: migration is enabled, but the required adjusted migration "
            << "schedules are missing or empty. Run src/param.py to generate both "
            << adjusted_standard_mig_json << " and "
            << adjusted_inverted_mig_json << ".\n";
        return EXIT_FAILURE;
    }
    if (migration_enabled &&
        (adjusted_standard_schedule.front().first > 0.0 ||
         adjusted_inverted_schedule.front().first > 0.0)) {
        std::cerr
            << "Error: adjusted migration schedules must contain a matrix active "
            << "at generation 0. Regenerate "
            << adjusted_standard_mig_json << " and "
            << adjusted_inverted_mig_json << " with src/param.py.\n";
        return EXIT_FAILURE;
    }

    // Horizontal SMC consumes the precomputed context-adjusted schedules.
    // Raw migration rates are retained for the initial x0 simulation, which
    // applies its existing q1/q2 conversion internally.
    const auto& horizontal_standard_schedule = adjusted_standard_schedule;
    const auto& horizontal_inverted_schedule = adjusted_inverted_schedule;

    vector<vector<double>> mig_prob;
    vector<vector<double>> mig_prob_cut;
    size_t next_idx = 0;

    auto resetMigrationState = [&]() {
        mig_prob.clear();
        mig_prob_cut.clear();
        next_idx = 0;
        if (schedule.empty()) {
            mig_prob = buildMigMatrix(params);
            mig_prob_cut = mig_prob;
            return;
        }

        double g0 = 0.0;
        while (next_idx < schedule.size() && schedule[next_idx].first <= g0) {
            mig_prob = adaptMatrixForPops(schedule[next_idx].second, params.paramData->popSizeVec.size());
            ++next_idx;
        }
        if (mig_prob.empty()) {
            mig_prob = adaptMatrixForPops(schedule.front().second, params.paramData->popSizeVec.size());
        }
        // For the cut lineage, use the adjusted migration schedule if available.
        if (adjusted_standard_schedule.empty()) {
            mig_prob_cut = mig_prob;
        } else {

            for (const auto& item : adjusted_standard_schedule) {
                if (item.first <= g0) {
                    mig_prob_cut = adaptMatrixForPops(
                        item.second, params.paramData->popSizeVec.size());
                } else {
                    break;
                }
            }
        }
    };

    const string hopTracePath = writeAllDiagnostics
        ? pathJoin(output_dir, "smc_hop_trace.csv")
        : "smc_hop_trace.csv";
    const string hopEventsPath = writeAllDiagnostics
        ? pathJoin(output_dir, "smc_hop_events.csv")
        : "smc_hop_events.csv";
    std::ofstream hopTrace(hopTracePath);
    if (hopTrace.is_open()) {
        hopTrace << "run,hop,current_x,raw_delta_x,used_delta_x,next_x,rho,Li_sum,Ls_sum,root_time\n";
    }
    std::ofstream hopEvents(hopEventsPath);
    if (hopEvents.is_open()) {
        hopEvents << "run,hop,current_x,event,event_time,raw_delta_x,used_delta_x,next_x\n";
    }
    std::ofstream treeSnapshotsCSV;
    if (params.paramData->csvSnapshots) {
        treeSnapshotsCSV.open("smc_tree_snapshots.csv");
    }
    if (treeSnapshotsCSV.is_open()) {
        treeSnapshotsCSV << "run,hop,x_start,x_end,node_id,parent_id,time,pop,inversion\n";
    }
    std::ofstream treeSnapshotsBin;
    if (params.paramData->binarySnapshots) {
        treeSnapshotsBin.open("smc_tree_snapshots.bin", std::ios::binary);
    }
    if (treeSnapshotsBin.is_open()) {
        treeSnapshotsBin.write("SMCTREE1", 8);
    }

    for (int timer = 0; timer < (int)nRuns; ++timer)
    {
        resetMigrationState();
        params.setPhi();
        params.setSNPs();
        params.setCarriers();

        World *world = new World(params.getpData());
        while (!world->simulationFinished())
        {
            if (!schedule.empty() && next_idx < schedule.size()) {
                double now   = world->nGenerations();
                double tnext = schedule[next_idx].first;
                if (now >= tnext) {
                    mig_prob = adaptMatrixForPops(schedule[next_idx].second, params.paramData->popSizeVec.size());
                    ++next_idx;
                }
            }
            world->simulateGeneration_SMC(mig_prob);
        }

        vector<shared_ptr<ARGNode>> allNodes = world->getARGVec();
        if (nSites == 0 || allNodes.empty()) {
            std::cerr << "No sites or ARG nodes available.\n";
            delete world;
            break;
        }

        const double startX = params.paramData->smcRange.L;
        const double argStartX = startX / params.paramData->BasesPerMorgan;
        SiteNode geneTree(argStartX, allNodes.back());
        geneTree.calcBranchLengths(0);
        geneTree.calcBranchLengths_informative(0);
        if (writeAllDiagnostics) {
            geneTree.writeCSV(pathJoin(output_dir, "genetree_first_site.csv"));
            geneTree.writeDOT(pathJoin(tree_dir, "genetree_first_site.dot"));
            allNodes.back()->writeDOT(pathJoin(tree_dir, "argtree.dot"));
        }

        const double r = 1.0e-8;
        TreeNode* activeTree = buildX0TreeFromARGPreserveUnary(argStartX, allNodes.back());
        if (!activeTree) {
            std::cerr << "Could not build unary-preserving x0 tree from ARG.\n";
            delete world;
            break;
        }
        if (writeAllDiagnostics) {
            writeTreeDOT(activeTree, pathJoin(tree_dir, "genetree_x0_arg_unary.dot"));
            writeCollapsedTreeDOT(activeTree, pathJoin(tree_dir, "genetree_x0_arg_unary_collapsed.dot"));
        }
        unsigned long nextNodeId = getMaxId(activeTree) + 1;
        double currentX = startX;
        appendTreeSnapshotCSV(treeSnapshotsCSV, activeTree, timer, 0, currentX, currentX);
        appendTreeSnapshotBinary(treeSnapshotsBin, activeTree, timer, 0, currentX, currentX);
        vector<GeneFluxEvent_SMC> geneFluxActive;
        vector<GeneFluxEvent_SMC> geneFluxLog;
        vector<EdgeWeight> last_standard_edges;
        vector<EdgeWeight> last_inverted_edges;
        vector<bool> targetEmitted(params.paramData->targetSNPs.size(), false);
        if (targetMode) {
            const double eps = 1e-15;
            vector<double> x0Targets;
            for (size_t i = 0; i < params.paramData->targetSNPs.size(); ++i) {
                if (targetEmitted[i]) continue;
                const double targetX = params.paramData->targetSNPs[i];
                if (std::abs(targetX - currentX) <= eps) {
                    x0Targets.push_back(targetX);
                    targetEmitted[i] = true;
                }
            }
            for (double targetX : x0Targets) {
                std::ostringstream targetBase;
                targetBase << "genetree_target"
                           << formatCoordForFilename(targetX)
                           << "_hop0_x" << formatCoordForFilename(currentX);
                writeTreeArtifacts(activeTree, targetBase.str());
            }
        }
        int hop = 0;
        while (currentX < params.paramData->smcRange.R) {
            vector<EdgeWeight> standard_edges;
            vector<EdgeWeight> inverted_edges;
            collectEdgeWeightsFromTree(activeTree,
                                       *params.paramData,
                                       r,
                                       standard_edges,
                                       inverted_edges);
            last_standard_edges = standard_edges;
            last_inverted_edges = inverted_edges;
            double Ls_sum = 0.0, Li_sum = 0.0;
            for (const auto &e : standard_edges) Ls_sum += e.weight;
            for (const auto &e : inverted_edges) Li_sum += e.weight;
            const double rho = Ls_sum + Li_sum;
            if (rho <= 0.0) {
                break;
            }
            std::unordered_map<unsigned long, TreeNode*> clonedNodes;
            if (!treeShapeHeaderWritten) {
                writeTreeShapeHeader(tree_shape_path);
                treeShapeHeaderWritten = true;
            }
            appendTreeShapeRow(tree_shape_path,
                               timer,
                               hop,
                               currentX,
                               rho,
                               summarizeTreeShape(activeTree));

            TreeNode* workingTree = cloneTreeWithMap(activeTree, clonedNodes);
            unsigned long cutParentId = 0;
            unsigned long cutChildId = 0;
            bool picked = pickWeightedEdge(standard_edges, inverted_edges, &cutParentId, &cutChildId);
            if (!picked) {
                freeTree(workingTree);
                std::cerr << "SMC cut-tree step failed (could not sample weighted edge).\n";
                break;
            }
            TreeNode* p = clonedNodes.count(cutParentId) ? clonedNodes[cutParentId] : nullptr;
            TreeNode* c = clonedNodes.count(cutChildId) ? clonedNodes[cutChildId] : nullptr;
            if (!p || !c || c->parent != p) {
                freeTree(workingTree);
                std::cerr << "SMC cut-tree step failed (sampled edge not found in tree).\n";
                break;
            }
            TreeNode* cutSubtree = nullptr;
            double cutStartTime = 0.0;
            bool cutOk = cutEdgeRandomWithCleanup(workingTree, p, c, &cutSubtree, &cutStartTime);
            if (!cutOk) {
                freeTree(workingTree);
                std::cerr << "SMC cut-tree step skipped (invalid cut edge).\n";
                break;
            }
            SMCStepOutcome outcome;
            bool ok = simulateSMCOnTree_SMC(workingTree,
                                            cutSubtree,
                                            cutStartTime,
                                            *params.paramData,
                                            mig_prob_cut,
                                            horizontal_standard_schedule,
                                            horizontal_inverted_schedule,
                                            currentX,
                                            hop,
                                            &outcome,
                                            timer,
                                            nextNodeId,
                                            hopEvents);

            // Append every vertical rate calc during this hop
            if (!outcome.coalescenceRows.empty()) {
                std::ofstream coalescenceDiagnostics(
                    coalescence_diagnostic_path, std::ios::app);
                if (coalescenceDiagnostics.is_open()) {
                    for (const auto& row : outcome.coalescenceRows) {
                        coalescenceDiagnostics << timer << "," << row;
                    }
                }
            }

            for (const auto& evt : outcome.geneFluxEvents) {
                geneFluxActive.push_back(evt);
            }

            if (!ok) {
                if (hopEvents.is_open()) {
                    for (const auto& row : outcome.eventRows) {
                        hopEvents << timer << "," << row;
                    }
                    hopEvents << timer << ","
                               << hop << ","
                               << currentX << ","
                               << "reattach_failed,"
                               << ",,,\n";
                }
                freeTree(workingTree);
                std::cerr << "SMC cut-tree step failed to reattach.\n";
                break;
            }

            trimUnaryRootStem(workingTree);
            freeTree(activeTree);
            activeTree = workingTree;

            const double rawHopDelta = sampleHopDeltaFromRho(rho);
            if (rawHopDelta <= 0.0) {
                break;
            }
            double hopDelta = rawHopDelta;
            if (currentX + hopDelta > params.paramData->smcRange.R) {
                hopDelta = params.paramData->smcRange.R - currentX;
            }
            if (hopDelta <= 0.0) {
                break;
            }
            double nextX = currentX + hopDelta;

            const double activeEps = 1e-12;
            geneFluxActive.erase(
                std::remove_if(geneFluxActive.begin(), geneFluxActive.end(),
                               [currentX, activeEps](const GeneFluxEvent_SMC& evt) {
                                   return evt.endX <= currentX + activeEps;
                               }),
                geneFluxActive.end());

            if (!geneFluxActive.empty()) {
                size_t minIdx = 0;
                for (size_t i = 1; i < geneFluxActive.size(); ++i) {
                    if (geneFluxActive[i].endX < geneFluxActive[minIdx].endX) {
                        minIdx = i;
                    }
                }
                const GeneFluxEvent_SMC minEvt = geneFluxActive[minIdx];
                if (nextX >= minEvt.endX) {
                    nextX = minEvt.endX;
                    geneFluxLog.push_back(minEvt);
                    geneFluxActive.erase(geneFluxActive.begin() + static_cast<long>(minIdx));
                }
            }
            const double finalHopDelta = nextX - currentX;
            if (finalHopDelta <= 0.0) {
                continue;
            }
            appendTreeSnapshotCSV(treeSnapshotsCSV, activeTree, timer, hop + 1, currentX, nextX);
            appendTreeSnapshotBinary(treeSnapshotsBin, activeTree, timer, hop + 1, currentX, nextX);

            bool writeThisHop = writeAllDiagnostics;
            vector<double> targetsForThisHop;
            if (targetMode) {
                const double eps = 1e-15;
                for (size_t i = 0; i < params.paramData->targetSNPs.size(); ++i) {
                    if (targetEmitted[i]) continue;
                    const double targetX = params.paramData->targetSNPs[i];
                    if (targetX < currentX - eps) {
                        targetEmitted[i] = true;
                        continue;
                    }
                    if (targetX <= nextX + eps) {
                        writeThisHop = true;
                        targetsForThisHop.push_back(targetX);
                        targetEmitted[i] = true;
                    }
                }
            }

            if (writeThisHop) {
                vector<std::string> artifactBases;
                if (targetMode) {
                    for (double targetX : targetsForThisHop) {
                        std::ostringstream targetBase;
                        targetBase << "genetree_target"
                                   << formatCoordForFilename(targetX)
                                   << "_hop" << (hop + 1)
                                   << "_x" << formatCoordForFilename(currentX)
                                   << "_to_" << formatCoordForFilename(nextX);
                        artifactBases.push_back(targetBase.str());
                    }
                } else {
                    std::ostringstream hopBase;
                    hopBase << "genetree_hop" << (hop + 1);
                    artifactBases.push_back(hopBase.str());
                }

                for (const auto& hopBase : artifactBases) {
                    const string hopBasePath = pathJoin(tree_dir, hopBase);
                    writeTreeArtifacts(activeTree, hopBasePath);

                    if (targetMode) {
                        writeEdgeWeightsCSV(standard_edges, hopBasePath + "_edge_weights_standard.csv");
                        writeEdgeWeightsCSV(inverted_edges, hopBasePath + "_edge_weights_inverted.csv");
                    }
                }
            }

            if (hopTrace.is_open()) {
                hopTrace << std::setprecision(17);
                hopTrace << timer << ","
                         << hop << ","
                         << currentX << ","
                         << rawHopDelta << ","
                         << finalHopDelta << ","
                         << nextX << ","
                         << rho << ","
                         << Li_sum << ","
                         << Ls_sum << ","
                         << (activeTree ? activeTree->time : 0.0) << "\n";
            }
            if (hopEvents.is_open()) {
                for (const auto& row : outcome.eventRows) {
                    hopEvents << timer << "," << row;
                }
                hopEvents << timer << ","
                          << hop << ","
                          << currentX << ","
                          << "hop_summary,"
                          << ","
                          << rawHopDelta << ","
                          << finalHopDelta << ","
                          << nextX << "\n";
            }
            currentX = nextX;
            ++hop;
        }

        if (writeAllDiagnostics) {
            std::ofstream gf_active(pathJoin(output_dir, "smc_gene_flux_active.csv"));
            if (gf_active.is_open()) {
                gf_active << "x_start,x_end,node_id,type\n";
                for (const auto& evt : geneFluxActive) {
                    gf_active << evt.startX << "," << evt.endX << "," << evt.nodeId
                              << "," << evt.type << "\n";
                }
            }
        }
        if (writeAllDiagnostics) {
            std::ofstream gf_log(pathJoin(output_dir, "smc_gene_flux_log.csv"));
            if (gf_log.is_open()) {
                gf_log << "x_start,x_end,node_id,type\n";
                for (const auto& evt : geneFluxLog) {
                    gf_log << evt.startX << "," << evt.endX << "," << evt.nodeId
                           << "," << evt.type << "\n";
                }
            }
        }

        hopEvents.flush();
        hopTrace.flush();

        if (writeAllDiagnostics) {
            writeTreeDOT(activeTree, pathJoin(tree_dir, "genetree_modified.dot"));
        }
        freeTree(activeTree);

        if (writeAllDiagnostics) {
            double std_len = 0.0;
            double inv_len = 0.0;
            geneTree.getTotalLengthByInversion(std_len, inv_len);
            std::cerr << "Ls(x) = " << std_len << "\n";
            std::cerr << "Li(x) = " << inv_len << "\n";
            std::cerr << "L(x) = " << (std_len + inv_len) << "\n";
            std::cerr << "L(x) (getTotalLength func) = "
                      << geneTree.getTotalLength(0.0) << "\n";
        }

        if (writeAllDiagnostics) {
            writeEdgeWeightsCSV(last_standard_edges, pathJoin(output_dir, "edge_weights_standard.csv"));
            writeEdgeWeightsCSV(last_inverted_edges, pathJoin(output_dir, "edge_weights_inverted.csv"));
        }

        delete world;
        printProgress(timer + 1, static_cast<int>(nRuns), progressBucket);
    }

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cerr << "Elapsed time: " << elapsed_seconds.count() << "s\n";

    return 0;
}
