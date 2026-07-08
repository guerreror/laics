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
    const vector<unsigned int>& pop_sizes,
    const vector<double>& inv_freqs,
    double r,
    vector<EdgeWeight>& standard_edges,
    vector<EdgeWeight>& inverted_edges)
{
    if (!node) return;
    for (auto* child : node->children) {
        if (!child) continue;
        const unsigned int pop = child->context.pop;
        if (pop < pop_sizes.size() && pop < inv_freqs.size()) {
            const double popN = pop_sizes[pop];
            const double pI = inv_freqs[pop];
            const double branchL = node->time - child->time;
            if (branchL > 0.0) {
                const double base = 2.0 * r * popN * branchL;
                if (child->context.inversion == 1) {
                    inverted_edges.push_back({node->id, child->id, base * pI});
                } else {
                    standard_edges.push_back({node->id, child->id, base * (1.0 - pI)});
                }
            }
        }
        collectEdgeWeightsFromTree(child, pop_sizes, inv_freqs, r, standard_edges, inverted_edges);
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

    vector<string> param_vec(argv + 2, argv + argc);
    Parameters params(infile.str().c_str(), param_vec);

    unsigned int nRuns = params.paramData->nRuns;
    unsigned int nSites = params.paramData->n_SNPs;
    const bool targetMode = !params.paramData->targetSNPs.empty();
    const bool writeAllDiagnostics = !targetMode && params.paramData->smcVerbose;

    const std::string mig_json = "src/migration_matrices.json";
    auto schedule = readMigrationSchedule(mig_json);

    vector<vector<double>> mig_prob;
    vector<vector<double>> mig_prob_cut;
    size_t next_idx = 0;

    if (!schedule.empty()) {
        double g0 = 0.0;
        while (next_idx < schedule.size() && schedule[next_idx].first <= g0) {
            mig_prob = adaptMatrixForPops(schedule[next_idx].second, params.paramData->popSizeVec.size());
            ++next_idx;
        }
        if (mig_prob.empty()) {
            mig_prob = adaptMatrixForPops(schedule.front().second, params.paramData->popSizeVec.size());
        }
        mig_prob_cut = mig_prob;
    } else {
        mig_prob = buildMigMatrix(params);
        mig_prob_cut = mig_prob;
    }

    for (int timer = 0; timer < (int)nRuns; ++timer)
    {
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
            geneTree.writeCSV("genetree_first_site.csv");
            geneTree.writeDOT("genetree_first_site.dot");
            allNodes.back()->writeDOT("argtree.dot");
        }

        const double r = 1.0e5;
        TreeNode* activeTree = buildX0TreeFromARGPreserveUnary(argStartX, allNodes.back());
        if (!activeTree) {
            std::cerr << "Could not build unary-preserving x0 tree from ARG.\n";
            delete world;
            break;
        }
        if (writeAllDiagnostics) {
            writeTreeDOT(activeTree, "genetree_x0_arg_unary.dot");
            writeCollapsedTreeDOT(activeTree, "genetree_x0_arg_unary_collapsed.dot");
        }
        double currentX = startX;
        vector<GeneFluxEvent_SMC> geneFluxActive;
        vector<GeneFluxEvent_SMC> geneFluxLog;
        vector<EdgeWeight> last_standard_edges;
        vector<EdgeWeight> last_inverted_edges;
        vector<bool> targetEmitted(params.paramData->targetSNPs.size(), false);
        bool hopTraceHeaderWritten = false;
        bool hopEventsHeaderWritten = false;
        if (writeAllDiagnostics) {
            std::ofstream hoplog("smc_hop_trace.csv");
            if (hoplog.is_open()) {
                hoplog << "hop,current_x,raw_delta_x,used_delta_x,next_x,rho,Li_sum,Ls_sum\n";
                hopTraceHeaderWritten = true;
            }
            std::ofstream hop_events("smc_hop_events.csv");
            if (hop_events.is_open()) {
                hop_events << "hop,current_x,event,event_time,raw_delta_x,used_delta_x,next_x\n";
                hopEventsHeaderWritten = true;
            }
        }

        auto ensureHopTraceHeader = [&]() {
            if (hopTraceHeaderWritten) return;
            std::ofstream hoplog("smc_hop_trace.csv");
            if (hoplog.is_open()) {
                hoplog << "hop,current_x,raw_delta_x,used_delta_x,next_x,rho,Li_sum,Ls_sum\n";
                hopTraceHeaderWritten = true;
            }
        };

        auto ensureHopEventsHeader = [&]() {
            if (hopEventsHeaderWritten) return;
            std::ofstream hop_events("smc_hop_events.csv");
            if (hop_events.is_open()) {
                hop_events << "hop,current_x,event,event_time,raw_delta_x,used_delta_x,next_x\n";
                hopEventsHeaderWritten = true;
            }
        };

        const int maxHopsSkeleton = 1000;
        for (int hop = 0; hop < maxHopsSkeleton; ++hop) {
            vector<EdgeWeight> standard_edges;
            vector<EdgeWeight> inverted_edges;
            collectEdgeWeightsFromTree(activeTree,
                                       params.paramData->popSizeVec,
                                       params.paramData->initialFreqs,
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

            TreeNode* workingTree = cloneTree(activeTree);
            unsigned long cutParentId = 0;
            unsigned long cutChildId = 0;
            bool picked = pickWeightedEdge(standard_edges, inverted_edges, &cutParentId, &cutChildId);
            if (!picked) {
                freeTree(workingTree);
                std::cerr << "SMC cut-tree step failed (could not sample weighted edge).\n";
                break;
            }
            TreeNode* p = findNodeById(workingTree, cutParentId);
            TreeNode* c = findNodeById(workingTree, cutChildId);
            if (!p || !c || c->parent != p) {
                freeTree(workingTree);
                std::cerr << "SMC cut-tree step failed (sampled edge not found in tree).\n";
                break;
            }
            TreeNode* cutSubtree = nullptr;
            unsigned long maxId = 0;
            std::function<void(TreeNode*)> gather = [&](TreeNode* n){
                if (!n) return;
                if (n->id > maxId) maxId = n->id;
                for (auto* ch : n->children) gather(ch);
            };
            gather(workingTree);
            const unsigned long cutpointId = maxId + 1;
            bool cutOk = cutEdgeRandomWithCleanup(workingTree, p, c, cutpointId, &cutSubtree);
            if (!cutOk) {
                freeTree(workingTree);
                std::cerr << "SMC cut-tree step skipped (invalid cut edge).\n";
                break;
            }

            SMCStepOutcome outcome;
            bool ok = simulateSMCOnTree_SMC(workingTree,
                                            cutSubtree,
                                            *params.paramData,
                                            mig_prob_cut,
                                            currentX,
                                            hop,
                                            &outcome,
                                            "");

            for (const auto& evt : outcome.geneFluxEvents) {
                geneFluxActive.push_back(evt);
            }

            if (!ok) {
                freeTree(workingTree);
                std::cerr << "SMC cut-tree step failed to reattach.\n";
                break;
            }

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

            bool writeThisHop = writeAllDiagnostics;
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
                        targetEmitted[i] = true;
                    }
                }
            }

            if (writeThisHop) {
                std::ostringstream hopBase;
                hopBase << "genetree_hop" << (hop + 1);
                writeTreeArtifacts(activeTree, hopBase.str());

                if (targetMode) {
                    writeEdgeWeightsCSV(standard_edges, hopBase.str() + "_edge_weights_standard.csv");
                    writeEdgeWeightsCSV(inverted_edges, hopBase.str() + "_edge_weights_inverted.csv");
                }
            }

            if (writeThisHop) {
                ensureHopTraceHeader();
                std::ofstream hoplog("smc_hop_trace.csv", std::ios::app);
                if (hoplog.is_open()) {
                    hoplog << hop << ","
                           << currentX << ","
                           << rawHopDelta << ","
                           << hopDelta << ","
                           << nextX << ","
                           << rho << ","
                           << Li_sum << ","
                           << Ls_sum << "\n";
                }
            }
            if (writeThisHop) {
                ensureHopEventsHeader();
                std::ofstream hop_events("smc_hop_events.csv", std::ios::app);
                if (hop_events.is_open()) {
                    for (const auto& row : outcome.eventRows) {
                        hop_events << row;
                    }
                    hop_events << hop << ","
                               << currentX << ","
                               << "hop_summary,"
                               << ","
                               << rawHopDelta << ","
                               << hopDelta << ","
                               << nextX << "\n";
                }
            }
            currentX = nextX;
        }

        if (writeAllDiagnostics) {
            std::ofstream gf_active("smc_gene_flux_active.csv");
            if (gf_active.is_open()) {
                gf_active << "x_start,x_end,node_id,type\n";
                for (const auto& evt : geneFluxActive) {
                    gf_active << evt.startX << "," << evt.endX << "," << evt.nodeId
                              << "," << evt.type << "\n";
                }
            }
        }
        if (writeAllDiagnostics) {
            std::ofstream gf_log("smc_gene_flux_log.csv");
            if (gf_log.is_open()) {
                gf_log << "x_start,x_end,node_id,type\n";
                for (const auto& evt : geneFluxLog) {
                    gf_log << evt.startX << "," << evt.endX << "," << evt.nodeId
                           << "," << evt.type << "\n";
                }
            }
        }

        if (writeAllDiagnostics) {
            writeTreeDOT(activeTree, "genetree_modified.dot");
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
            writeEdgeWeightsCSV(last_standard_edges, "edge_weights_standard.csv");
            writeEdgeWeightsCSV(last_inverted_edges, "edge_weights_inverted.csv");
        }

        delete world;
        break;
    }

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cerr << "Elapsed time: " << elapsed_seconds.count() << "s\n";

    return 0;
}
