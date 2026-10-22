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
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <future>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <sys/stat.h>

#include <tskit/tables.h>

using namespace std;

#include "parameters.h"
#include "world.h"
#include "ran_mk.h"
#include "migprob.h"
#include "sitenode.h"
#include "treemod.h"
#include "simulate_smc.h"
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
                             TreeNode** outParentNode,
                             TreeNode** outChildNode)
{
    if (!outParentNode || !outChildNode) return false;
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
            *outParentNode = e.parentNode;
            *outChildNode = e.childNode;
            return true;
        }
    }
    *outParentNode = all.back().parentNode;
    *outChildNode = all.back().childNode;
    return true;
}

struct GeneFluxBoundaryRef {
    double endX;
    unsigned long nodeId;
};

struct GeneFluxBoundaryEarlier {
    bool operator()(const GeneFluxBoundaryRef& a, const GeneFluxBoundaryRef& b) const {
        return a.endX > b.endX;
    }
};

struct ActiveGeneFluxSet {
    std::unordered_map<unsigned long, GeneFluxEvent_SMC> byId;
    std::priority_queue<GeneFluxBoundaryRef,
                        std::vector<GeneFluxBoundaryRef>,
                        GeneFluxBoundaryEarlier> byEnd;

    void prune() {
        while (!byEnd.empty()) {
            const auto& top = byEnd.top();
            auto it = byId.find(top.nodeId);
            if (it != byId.end() && std::abs(it->second.endX - top.endX) <= 1e-9) break;
            byEnd.pop();
        }
    }

    size_t size() const { return byId.size(); }

    void addIncoming(const vector<GeneFluxEvent_SMC>& incoming, double currentX, double eps) {
        for (const auto& evt : incoming) {
            if (evt.endX <= currentX + eps) continue;
            byId[evt.nodeId] = evt;
            byEnd.push({evt.endX, evt.nodeId});
        }
    }

    bool firstExpired(double currentX, double eps, GeneFluxEvent_SMC& out) {
        prune();
        if (byEnd.empty() || byEnd.top().endX > currentX + eps) return false;
        out = byId[byEnd.top().nodeId];
        return true;
    }

    bool nextEnd(double& endX) {
        prune();
        if (byEnd.empty()) return false;
        endX = byEnd.top().endX;
        return true;
    }

    void erase(unsigned long nodeId) {
        byId.erase(nodeId);
    }

    void removeDeleted(const vector<unsigned long>& deletedIds) {
        for (unsigned long id : deletedIds) erase(id);
    }
};

static bool processGeneFluxBoundariesAtX(
    TreeNode*& activeTree,
    ActiveGeneFluxSet& geneFluxActive,
    const Parameters::ParameterData& paramData,
    const vector<vector<double>>& mig_prob_cut,
    const vector<pair<double, vector<vector<double>>>>& standard_mig_schedule,
    const vector<pair<double, vector<vector<double>>>>& inverted_mig_schedule,
    double currentX,
    int hop,
    int timer,
    unsigned long& nextNodeId,
    std::ofstream& hopEvents)
{
    const double eps = 1e-12;
    while (true) {
        GeneFluxEvent_SMC evt;
        if (!geneFluxActive.firstExpired(currentX, eps, evt)) return true;

        TreeNode* gfNode = findNodeById(activeTree, evt.nodeId);
        if (!gfNode) {
            geneFluxActive.erase(evt.nodeId);
            if (hopEvents.is_open()) {
                hopEvents << timer << "," << hop << "," << currentX
                          << ",gene_flux_boundary_stale,,,,\n";
            }
            continue;
        }
        if (!gfNode->parent || gfNode->children.size() != 1) {
            if (hopEvents.is_open()) {
                hopEvents << timer << "," << hop << "," << currentX
                          << ",gene_flux_boundary_invalid," << gfNode->time << ",,,\n";
            }
            return false;
        }

        TreeNode* cutRoot = nullptr;
        vector<unsigned long> deletedIds;
        const double restartTime = gfNode->time;
        TreeNode* cutNode = gfNode->children.front();
        const bool cutOk = cutAtNodeWithUnaryCleanupCollect(activeTree,
                                                            cutNode,
                                                            restartTime,
                                                            &cutRoot,
                                                            deletedIds);
        if (!cutOk || !activeTree || !cutRoot) {
            if (hopEvents.is_open()) {
                hopEvents << timer << "," << hop << "," << currentX
                          << ",gene_flux_boundary_cut_failed," << restartTime << ",,,\n";
            }
            return false;
        }

        geneFluxActive.removeDeleted(deletedIds);
        if (hopEvents.is_open()) {
            hopEvents << timer << "," << hop << "," << currentX
                      << ",gene_flux_boundary," << restartTime << ",,,\n";
        }

        SMCStepOutcome outcome;
        const bool ok = simulateSMCOnTree_SMC(activeTree,
                                              cutRoot,
                                              restartTime,
                                              paramData,
                                              mig_prob_cut,
                                              standard_mig_schedule,
                                              inverted_mig_schedule,
                                              currentX,
                                              hop,
                                              &outcome,
                                              timer,
                                              nextNodeId,
                                              hopEvents);
        geneFluxActive.addIncoming(outcome.geneFluxEvents, currentX, eps);
        if (!ok) {
            if (hopEvents.is_open()) {
                for (const auto& row : outcome.eventRows) {
                    hopEvents << timer << "," << row;
                }
                hopEvents << timer << "," << hop << "," << currentX
                          << ",gene_flux_boundary_reattach_failed," << restartTime << ",,,\n";
            }
            return false;
        }
        trimUnaryRootStem(activeTree);
    }
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
        const SMCActiveState& state = activeStateAtTime_SMC(params, child->time);
        if (leavesRetainedTree && pop < state.popSizes.size() && pop < state.invFreqs.size()) {
            const double pI = state.invFreqs[pop];
            const double branchL = node->time - child->time;
            if (branchL > 0.0) {
                const double base = r * branchL;
                if (child->context.inversion == 1) {
                    inverted_edges.push_back({node->id, child->id, base * pI, node, child});
                } else {
                    standard_edges.push_back({node->id, child->id, base * (1.0 - pI), node, child});
                }
            }
        }
        collectEdgeWeightsFromTree(child, params, r, standard_edges, inverted_edges);
    }
}

struct TskitCNodeRecord {
    int32_t run;
    uint64_t id;
    double time;
    uint32_t pop;
};

struct TskitCEdgeInterval {
    int32_t run;
    double left;
    double right;
    uint64_t parent;
    uint64_t child;
};

struct TskitCEdgeKey {
    uint64_t parent;
    uint64_t child;
    bool operator==(const TskitCEdgeKey& other) const {
        return parent == other.parent && child == other.child;
    }
};

struct TskitCKeyHash {
    size_t operator()(const TskitCEdgeKey& k) const {
        uint64_t h = k.parent;
        h ^= k.child + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return std::hash<uint64_t>()(h);
    }
};

struct TskitCOpenEdge {
    size_t index;
    uint64_t generation;
};

struct TskitCStore {
    std::vector<TskitCNodeRecord> nodes;
    std::vector<TskitCEdgeInterval> edges;
    std::vector<uint8_t> seenNodeIds;
    std::unordered_map<TskitCEdgeKey, TskitCOpenEdge, TskitCKeyHash> openEdges;
    uint64_t edgeGeneration = 0;

    static void checkTsk(int ret, const char* where) {
        if (ret < 0) {
            throw std::runtime_error(std::string(where) + ": " + tsk_strerror(ret));
        }
    }

    void collect(TreeNode* node, int run, double left, double right) {
        if (!node) return;
        const uint64_t nodeId = static_cast<uint64_t>(node->id);
        if (nodeId >= seenNodeIds.size()) seenNodeIds.resize(nodeId + 1, 0);
        if (!seenNodeIds[nodeId]) {
            seenNodeIds[nodeId] = 1;
            nodes.push_back({static_cast<int32_t>(run), nodeId, node->time,
                             static_cast<uint32_t>(node->context.pop)});
        }
        for (auto* child : node->children) {
            if (!child) continue;
            const uint64_t childId = static_cast<uint64_t>(child->id);
            const TskitCEdgeKey ek{nodeId, childId};
            if (right > left) {
                auto it = openEdges.find(ek);
                if (it != openEdges.end() &&
                    std::abs(edges[it->second.index].right - left) <= 1e-9) {
                    edges[it->second.index].right = right;
                } else {
                    edges.push_back({static_cast<int32_t>(run), left, right, nodeId, childId});
                    if (it == openEdges.end()) {
                        it = openEdges.emplace(ek, TskitCOpenEdge{}).first;
                    }
                    it->second.index = edges.size() - 1;
                }
                it->second.generation = edgeGeneration;
            }
            collect(child, run, left, right);
        }
    }

    void append(TreeNode* tree, int run, double left, double right) {
        if (!tree) return;
        ++edgeGeneration;
        collect(tree, run, left, right);
        for (auto it = openEdges.begin(); it != openEdges.end();) {
            if (it->second.generation != edgeGeneration) {
                it = openEdges.erase(it);
            } else {
                ++it;
            }
        }
    }

    void finishRun(int run) {
        (void) run;
        openEdges.clear();
        seenNodeIds.clear();
        edgeGeneration = 0;
    }

    size_t storedBytes() const {
        return nodes.size() * sizeof(TskitCNodeRecord)
             + edges.size() * sizeof(TskitCEdgeInterval);
    }

    void clear() {
        nodes.clear();
        edges.clear();
        finishRun(0);
    }

    void writeRun(const std::vector<const TskitCNodeRecord*>& runNodes,
                  const std::vector<const TskitCEdgeInterval*>& runEdges,
                  const std::string& outPath) const {
        uint64_t maxNodeId = 0;
        for (const auto* n : runNodes) maxNodeId = std::max(maxNodeId, n->id);
        std::vector<const TskitCNodeRecord*> allNodes(maxNodeId + 1, nullptr);
        for (const auto* n : runNodes) allNodes[n->id] = n;
        std::vector<TskitCNodeRecord> usedNodes(maxNodeId + 1);
        std::vector<uint8_t> used(maxNodeId + 1, 0);
        double sequenceLength = 0.0;
        uint32_t maxPop = 0;
        for (const auto* edge : runEdges) {
            const auto& e = *edge;
            if (e.right <= e.left) continue;
            sequenceLength = std::max(sequenceLength, e.right);
            if (e.parent < allNodes.size() && allNodes[e.parent]) {
                usedNodes[e.parent] = *allNodes[e.parent];
                used[e.parent] = 1;
                maxPop = std::max(maxPop, allNodes[e.parent]->pop);
            }
            if (e.child < allNodes.size() && allNodes[e.child]) {
                usedNodes[e.child] = *allNodes[e.child];
                used[e.child] = 1;
                maxPop = std::max(maxPop, allNodes[e.child]->pop);
            }
        }
        if (sequenceLength <= 0.0 || runNodes.empty()) return;
        std::vector<std::vector<uint64_t>> parentsByChild(maxNodeId + 1);
        std::vector<size_t> remainingChildren(maxNodeId + 1, 0);
        for (const auto* edge : runEdges) {
            const auto& e = *edge;
            if (e.right <= e.left) continue;
            parentsByChild[e.child].push_back(e.parent);
            ++remainingChildren[e.parent];
        }
        std::queue<uint64_t> ready;
        for (uint64_t id = 0; id <= maxNodeId; ++id) {
            if (used[id] && remainingChildren[id] == 0) ready.push(id);
        }
        while (!ready.empty()) {
            const uint64_t child = ready.front();
            ready.pop();
            for (const uint64_t parent : parentsByChild[child]) {
                if (usedNodes[parent].time <= usedNodes[child].time) {
                    usedNodes[parent].time = usedNodes[child].time + 1e-9;
                }
                if (--remainingChildren[parent] == 0) ready.push(parent);
            }
        }

        tsk_table_collection_t tables;
        checkTsk(tsk_table_collection_init(&tables, 0), "tsk_table_collection_init");
        tables.sequence_length = sequenceLength;
        const char timeUnits[] = "generations";
        checkTsk(tsk_table_collection_set_time_units(
                     &tables, timeUnits, static_cast<tsk_size_t>(strlen(timeUnits))),
                 "tsk_table_collection_set_time_units");

        for (uint32_t pop = 0; pop <= maxPop; ++pop) {
            checkTsk(tsk_population_table_add_row(&tables.populations, NULL, 0),
                     "tsk_population_table_add_row");
        }

        const size_t usedCount = static_cast<size_t>(
            std::count(used.begin(), used.end(), static_cast<uint8_t>(1)));
        std::vector<tsk_flags_t> nodeFlags;
        std::vector<double> nodeTimes;
        std::vector<tsk_id_t> nodePops;
        nodeFlags.reserve(usedCount);
        nodeTimes.reserve(usedCount);
        nodePops.reserve(usedCount);
        std::vector<tsk_id_t> rowFor(maxNodeId + 1, TSK_NULL);
        for (uint64_t id = 0; id <= maxNodeId; ++id) {
            if (!used[id]) continue;
            const auto& n = usedNodes[id];
            rowFor[id] = static_cast<tsk_id_t>(nodeTimes.size());
            nodeFlags.push_back((std::abs(n.time) <= 1e-12) ? TSK_NODE_IS_SAMPLE : 0);
            nodeTimes.push_back(n.time);
            nodePops.push_back(static_cast<tsk_id_t>(n.pop));
        }
        checkTsk(tsk_node_table_set_columns(&tables.nodes,
                     static_cast<tsk_size_t>(nodeTimes.size()), nodeFlags.data(),
                     nodeTimes.data(), nodePops.data(), NULL, NULL, NULL),
                 "tsk_node_table_set_columns");

        std::vector<double> edgeLeft;
        std::vector<double> edgeRight;
        std::vector<tsk_id_t> edgeParent;
        std::vector<tsk_id_t> edgeChild;
        edgeLeft.reserve(runEdges.size());
        edgeRight.reserve(runEdges.size());
        edgeParent.reserve(runEdges.size());
        edgeChild.reserve(runEdges.size());
        for (const auto* edge : runEdges) {
            const auto& e = *edge;
            if (e.right <= e.left) continue;
            edgeLeft.push_back(e.left);
            edgeRight.push_back(e.right);
            edgeParent.push_back(rowFor[e.parent]);
            edgeChild.push_back(rowFor[e.child]);
        }
        checkTsk(tsk_edge_table_set_columns(&tables.edges,
                     static_cast<tsk_size_t>(edgeLeft.size()), edgeLeft.data(), edgeRight.data(),
                     edgeParent.data(), edgeChild.data(), NULL, NULL),
                 "tsk_edge_table_set_columns");

        const char timestamp[] = "generated-by-laics-smc";
        const char record[] = "{\"software\":\"laics-smc\",\"writer\":\"tskit-c-api\"}";
        checkTsk(tsk_provenance_table_add_row(
                     &tables.provenances, timestamp, static_cast<tsk_size_t>(strlen(timestamp)),
                     record, static_cast<tsk_size_t>(strlen(record))),
                 "tsk_provenance_table_add_row");
        checkTsk(tsk_table_collection_sort(&tables, NULL, 0), "tsk_table_collection_sort");
        checkTsk(tsk_table_collection_build_index(&tables, 0), "tsk_table_collection_build_index");
        checkTsk(tsk_table_collection_dump(&tables, outPath.c_str(), 0), "tsk_table_collection_dump");
        tsk_table_collection_free(&tables);
    }

    void write(const std::string& outDir) const {
        mkdir(outDir.c_str(), 0755);
        std::unordered_map<int, std::vector<const TskitCNodeRecord*>> nodesByRun;
        std::unordered_map<int, std::vector<const TskitCEdgeInterval*>> edgesByRun;
        for (const auto& n : nodes) nodesByRun[n.run].push_back(&n);
        for (const auto& e : edges) edgesByRun[e.run].push_back(&e);

        std::vector<int> runs;
        runs.reserve(edgesByRun.size());
        for (const auto& kv : edgesByRun) runs.push_back(kv.first);
        std::sort(runs.begin(), runs.end());
        const unsigned int hw = std::thread::hardware_concurrency();
        const size_t maxWorkers = std::max<size_t>(1, std::min<size_t>(8, hw ? hw : 2));
        std::vector<std::future<void>> jobs;
        for (int run : runs) {
            jobs.push_back(std::async(std::launch::async, [this, &nodesByRun, &edgesByRun, outDir, run]() {
                writeRun(nodesByRun[run],
                         edgesByRun[run],
                         outDir + "/smc_tree_sequence_run" + std::to_string(run) + ".trees");
            }));
            if (jobs.size() >= maxWorkers) {
                jobs.front().get();
                jobs.erase(jobs.begin());
            }
        }
        for (auto& job : jobs) job.get();
    }
};

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
    std::ofstream details;
    std::streambuf* terminalErr = nullptr;
    if (const char* detailsPath = std::getenv("LAICS_RUN_DETAILS")) {
        details.open(detailsPath, std::ios::app);
        if (details) terminalErr = std::cerr.rdbuf(details.rdbuf());
    }

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
    if (terminalErr) {
        std::cerr.flush();
        std::cerr.rdbuf(terminalErr);
    }

    unsigned int nRuns = params.paramData->nRuns;
    unsigned int nSites = params.paramData->n_SNPs;
    const bool smcProfile = std::getenv("SMC_PROFILE") != nullptr;
    using ProfileClock = std::chrono::system_clock;
    auto elapsedSec = [](const ProfileClock::time_point& a,
                         const ProfileClock::time_point& b) {
        return std::chrono::duration<double>(b - a).count();
    };
    auto profNow = [&]() { return smcProfile ? ProfileClock::now() : ProfileClock::time_point{}; };
    auto profAdd = [&](double& acc, const ProfileClock::time_point& t0) {
        if (smcProfile) acc += elapsedSec(t0, ProfileClock::now());
    };
    double profArg = 0.0, profBuildX0 = 0.0, profBoundary = 0.0, profWeights = 0.0;
    double profCut = 0.0, profSim = 0.0, profAppend = 0.0, profTskWrite = 0.0;
    double profEdgePick = 0.0, profEdgeCopy = 0.0, profFluxBook = 0.0, profTreeReplace = 0.0;
    double profDelta = 0.0, profDiagnostics = 0.0;
    unsigned long long profActiveMinScans = 0, profFluxRemoveScans = 0, profMaxActiveFlux = 0;
    unsigned long long profHops = 0;

    int progressBucket = -1;
    printProgress(0, static_cast<int>(nRuns), progressBucket);

    const char* migPath = std::getenv("LAICS_MIGRATION_JSON");
    const std::string mig_json = migPath ? migPath : "migration_matrices.json";
    const char* adjustedStandardPath =
        std::getenv("LAICS_ADJUSTED_STANDARD_MIGRATION_JSON");
    const char* adjustedInvertedPath =
        std::getenv("LAICS_ADJUSTED_INVERTED_MIGRATION_JSON");
    const std::string adjusted_standard_mig_json = adjustedStandardPath
        ? adjustedStandardPath : "adjusted_standard_migration_matrices.json";
    const std::string adjusted_inverted_mig_json = adjustedInvertedPath
        ? adjustedInvertedPath : "adjusted_inverted_migration_matrices.json";
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
    const vector<pair<double, Matrix>> no_migration_schedule;
    const auto& horizontal_standard_schedule = migration_enabled
        ? adjusted_standard_schedule
        : no_migration_schedule;
    const auto& horizontal_inverted_schedule = migration_enabled
        ? adjusted_inverted_schedule
        : no_migration_schedule;

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
        // For the cut lineage, use the adjusted migration schedule only when
        // migration is enabled for this run.
        if (!migration_enabled) {
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

    std::ofstream hopTrace;
    if (params.paramData->hopInformation) {
        hopTrace.open("smc_hop_trace.csv");
    }
    if (hopTrace.is_open()) {
        hopTrace << "run,hop,current_x,raw_delta_x,used_delta_x,next_x,rho,Li_sum,Ls_sum,root_time\n";
    }
    std::ofstream hopEvents;
    if (params.paramData->hopInformation) {
        hopEvents.open("smc_hop_events.csv");
    }
    if (hopEvents.is_open()) {
        hopEvents << "run,hop,current_x,event,event_time,raw_delta_x,used_delta_x,next_x\n";
    }
    TskitCStore treeSequences;
    const size_t treeBatchBytes = size_t{512} * 1024 * 1024;
    const char* treeOutputPath = std::getenv("LAICS_TREE_OUTPUT_DIR");
    const std::string treeOutputDir =
        treeOutputPath && *treeOutputPath ? treeOutputPath : "smc_trees";

    for (int timer = 0; timer < (int)nRuns; ++timer)
    {
        resetMigrationState();
        params.setSNPs();
        params.setCarriers();

        auto profT0 = profNow();
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
        profAdd(profArg, profT0);

        vector<shared_ptr<ARGNode>> allNodes = world->getARGVec();
        if (nSites == 0 || allNodes.empty()) {
            std::cerr << "No sites or ARG nodes available.\n";
            delete world;
            break;
        }

        const double startX = params.paramData->smcRange.L;
        const double argStartX = startX / params.paramData->BasesPerMorgan;

        const double r = 1.0e-8;
        profT0 = profNow();
        TreeNode* activeTree = buildX0TreeFromARGPreserveUnary(argStartX, allNodes.back());
        profAdd(profBuildX0, profT0);
        if (!activeTree) {
            std::cerr << "Could not build unary-preserving x0 tree from ARG.\n";
            delete world;
            break;
        }
        unsigned long nextNodeId = getMaxId(activeTree) + 1;
        double currentX = startX;
        treeSequences.append(activeTree, timer, currentX, currentX);
        ActiveGeneFluxSet geneFluxActive;
        vector<EdgeWeight> standard_edges;
        vector<EdgeWeight> inverted_edges;

        int hop = 0;
        while (currentX < params.paramData->smcRange.R) {
            const double activeEps = 1e-12;
            profT0 = profNow();
            if (!processGeneFluxBoundariesAtX(activeTree,
                                             geneFluxActive,
                                             *params.paramData,
                                             mig_prob_cut,
                                             horizontal_standard_schedule,
                                             horizontal_inverted_schedule,
                                             currentX,
                                             hop,
                                             timer,
                                             nextNodeId,
                                             hopEvents)) {
                std::cerr << "SMC gene-flux boundary processing failed.\n";
                break;
            }
            profAdd(profBoundary, profT0);

            standard_edges.clear();
            inverted_edges.clear();
            profT0 = profNow();
            collectEdgeWeightsFromTree(activeTree,
                                       *params.paramData,
                                       r,
                                       standard_edges,
                                       inverted_edges);
            profAdd(profWeights, profT0);
            profT0 = profNow();
            double Ls_sum = 0.0, Li_sum = 0.0;
            for (const auto &e : standard_edges) Ls_sum += e.weight;
            for (const auto &e : inverted_edges) Li_sum += e.weight;
            profAdd(profEdgeCopy, profT0);
            const double rho = Ls_sum + Li_sum;
            if (rho <= 0.0) {
                break;
            }
            TreeNode* p = nullptr;
            TreeNode* c = nullptr;
            profT0 = profNow();
            bool picked = pickWeightedEdge(standard_edges,
                                           inverted_edges,
                                           &p,
                                           &c);
            profAdd(profEdgePick, profT0);
            if (!picked) {
                std::cerr << "SMC cut-tree step failed (could not sample weighted edge).\n";
                break;
            }
            if (!p || !c || c->parent != p) {
                std::cerr << "SMC cut-tree step failed (sampled edge not found in tree).\n";
                break;
            }
            TreeNode* workingTree = activeTree;
            TreeNode* cutSubtree = nullptr;
            double cutStartTime = 0.0;
            vector<unsigned long> deletedIds;
            profT0 = profNow();
            bool cutOk = cutEdgeRandomWithCleanupCollect(workingTree,
                                                         p,
                                                         c,
                                                         &cutSubtree,
                                                         &cutStartTime,
                                                         deletedIds);
            profAdd(profCut, profT0);
            if (!cutOk) {
                activeTree = workingTree;
                if (cutSubtree && !cutSubtree->parent) freeTree(cutSubtree);
                std::cerr << "SMC cut-tree step skipped (invalid cut edge).\n";
                break;
            }
            profT0 = profNow();
            profFluxRemoveScans += deletedIds.size();
            geneFluxActive.removeDeleted(deletedIds);
            profAdd(profFluxBook, profT0);
            SMCStepOutcome outcome;
            profT0 = profNow();
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
            profAdd(profSim, profT0);

            profT0 = profNow();
            geneFluxActive.addIncoming(outcome.geneFluxEvents, currentX, activeEps);
            if (geneFluxActive.size() > profMaxActiveFlux) profMaxActiveFlux = geneFluxActive.size();
            profAdd(profFluxBook, profT0);

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
                activeTree = workingTree;
                if (cutSubtree && !cutSubtree->parent) freeTree(cutSubtree);
                std::cerr << "SMC cut-tree step failed to reattach.\n";
                break;
            }

            profT0 = profNow();
            trimUnaryRootStem(workingTree);
            activeTree = workingTree;
            profAdd(profTreeReplace, profT0);

            profT0 = profNow();
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

            double nextFluxEnd = 0.0;
            if (geneFluxActive.nextEnd(nextFluxEnd)) {
                ++profActiveMinScans;
                if (nextX >= nextFluxEnd) {
                    nextX = nextFluxEnd;
                }
            }
            const double finalHopDelta = nextX - currentX;
            if (finalHopDelta <= 0.0) {
                profAdd(profDelta, profT0);
                continue;
            }
            profAdd(profDelta, profT0);
            profT0 = profNow();
            treeSequences.append(activeTree, timer, currentX, nextX);
            profAdd(profAppend, profT0);
            ++profHops;

            profT0 = profNow();
            {
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
            }
            {
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
            }
            currentX = nextX;
            ++hop;
            profAdd(profDiagnostics, profT0);
        }

        if (hopEvents.is_open()) hopEvents.flush();
        if (hopTrace.is_open()) hopTrace.flush();
        freeTree(activeTree);

        treeSequences.finishRun(timer);
        delete world;
        printProgress(timer + 1, static_cast<int>(nRuns), progressBucket);

        if (treeSequences.storedBytes() >= treeBatchBytes || timer + 1 == static_cast<int>(nRuns)) {
            const auto writeStart = ProfileClock::now();
            treeSequences.write(treeOutputDir);
            const double seconds = elapsedSec(writeStart, ProfileClock::now());
            if (smcProfile) profTskWrite += seconds;
            treeSequences.clear();
        }
    }

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cerr << "Total elapsed time: " << elapsed_seconds.count() << "s\n";
    if (smcProfile) {
        std::cerr << "SMC_PROFILE hops=" << profHops << "\n";
        std::cerr << "SMC_PROFILE arg_generation=" << profArg << "\n";
        std::cerr << "SMC_PROFILE build_x0=" << profBuildX0 << "\n";
        std::cerr << "SMC_PROFILE gene_flux_boundary=" << profBoundary << "\n";
        std::cerr << "SMC_PROFILE collect_edge_weights=" << profWeights << "\n";
        std::cerr << "SMC_PROFILE cut_cleanup=" << profCut << "\n";
        std::cerr << "SMC_PROFILE vertical_simulation=" << profSim << "\n";
        std::cerr << "SMC_PROFILE tsk_collect=" << profAppend << "\n";
        std::cerr << "SMC_PROFILE tsk_write=" << profTskWrite << "\n";
        std::cerr << "SMC_PROFILE edge_pick=" << profEdgePick << "\n";
        std::cerr << "SMC_PROFILE edge_copy_sum=" << profEdgeCopy << "\n";
        std::cerr << "SMC_PROFILE gene_flux_bookkeeping=" << profFluxBook << "\n";
        std::cerr << "SMC_PROFILE tree_replace_free=" << profTreeReplace << "\n";
        std::cerr << "SMC_PROFILE delta_and_boundaries=" << profDelta << "\n";
        std::cerr << "SMC_PROFILE diagnostics_checks=" << profDiagnostics << "\n";
        std::cerr << "SMC_PROFILE active_min_scan_items=" << profActiveMinScans << "\n";
        std::cerr << "SMC_PROFILE flux_remove_scan_items=" << profFluxRemoveScans << "\n";
        std::cerr << "SMC_PROFILE max_active_gene_flux=" << profMaxActiveFlux << "\n";
    }

    return 0;
}
