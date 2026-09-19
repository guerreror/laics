#include "smc_helpers.h"

#include <fstream>
#include <set>
#include <vector>
#include <functional>
#include <algorithm>

const SMCActiveState& activeStateAtTime_SMC(const Parameters::ParameterData& params, double t) {
    struct Event { double time; int type; size_t index; };
    static const Parameters::ParameterData* cachedParams = nullptr;
    static std::vector<double> times;
    static std::vector<SMCActiveState> states;

    if (cachedParams != &params) {
        cachedParams = &params;
        times.clear();
        states.clear();

        std::vector<Event> events;
    const size_t nPops = params.popSizeVec.size();
    events.reserve(params.demography.size() + params.speciation.size());

    if (!params.demography.empty() && params.demography[0] == 1) {
        const size_t stride = 1 + nPops;
        for (size_t i = 1; i + nPops < params.demography.size(); i += stride) {
            events.push_back({params.demography[i], 2, i});
        }
    }
    if (!params.speciation.empty() && params.speciation[0] == 1) {
        for (size_t i = 1; i + 4 < params.speciation.size(); i += 5) {
            events.push_back({params.speciation[i + 2], 1, i});
        }
    }

    std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
        if (a.time != b.time) return a.time < b.time;
        return a.type > b.type; // demography before speciation at the same time.
    });

    SMCActiveState state{params.popSizeVec, params.initialFreqs};
        times.push_back(-1.0);
        states.push_back(state);
    for (const auto& ev : events) {
        if (ev.type == 2) {
            for (size_t k = 0; k < state.popSizes.size() && ev.index + 1 + k < params.demography.size(); ++k) {
                const double size = params.demography[ev.index + 1 + k];
                if (size != 0) state.popSizes[k] = static_cast<unsigned int>(size);
            }
        } else {
            const unsigned int A = static_cast<unsigned int>(params.speciation[ev.index]);
            const unsigned int B = static_cast<unsigned int>(params.speciation[ev.index + 1]);
            if (A < state.popSizes.size() && B < state.popSizes.size()) {
                state.popSizes[A] = static_cast<unsigned int>(params.speciation[ev.index + 4]);
                state.invFreqs[A] = params.speciation[ev.index + 3];
                state.popSizes.erase(state.popSizes.begin() + static_cast<long>(B));
                state.invFreqs.erase(state.invFreqs.begin() + static_cast<long>(B));
            }
        }
            times.push_back(ev.time);
            states.push_back(state);
    }
    }

    const double query = t + 1e-9;
    const auto it = std::upper_bound(times.begin(), times.end(), query);
    const size_t idx = (it == times.begin()) ? 0 : static_cast<size_t>((it - times.begin()) - 1);
    return states[idx];
}

static bool segmentCarriesSite(double x, const std::vector<Segment>& segments) {
    for (const auto& s : segments) {
        if ((s.L < x && x < s.R) || s.L == x) {
            return true;
        }
    }
    return false;
}

static TreeNode* buildX0TreeRec(double x0, const std::shared_ptr<ARGNode>& argNode, TreeNode* parent) {
    if (!argNode || argNode->getTime() < 0) return nullptr;

    TreeNode* node = new TreeNode();
    node->id = argNode->getNodeNumber();
    node->time = argNode->getTime();
    node->context = argNode->getContext();
    node->parent = parent;

    const unsigned long nDesc = argNode->getNDescendants();
    for (unsigned long i = 0; i < nDesc; ++i) {
        const std::vector<Segment> segments = argNode->getDescendantSegmentVector(i);
        if (!segmentCarriesSite(x0, segments)) {
            continue;
        }

        TreeNode* child = buildX0TreeRec(x0, argNode->getDescendantNode(i), node);
        if (child) {
            node->children.push_back(child);
        }
    }

    return node;
}

TreeNode* buildX0TreeFromARGPreserveUnary(double x0, std::shared_ptr<ARGNode> argRoot) {
    return buildX0TreeRec(x0, argRoot, nullptr);
}

static bool isLeaf(const TreeNode* n) {
    return n && n->children.empty();
}

static bool isBranching(const TreeNode* n) {
    return n && n->children.size() > 1;
}

static bool shouldKeepNode(const TreeNode* n, const TreeNode* root) {
    if (!n) return false;
    if (n == root) return true;
    if (isLeaf(n)) return true;
    if (isBranching(n)) return true;
    return false; // unary internal node gets contracted
}

static void collectKeptNodes(TreeNode* root, std::set<TreeNode*>& kept) {
    if (!root) return;
    std::function<void(TreeNode*)> dfs = [&](TreeNode* n) {
        if (!n) return;
        if (shouldKeepNode(n, root)) {
            kept.insert(n);
        }
        for (auto* ch : n->children) dfs(ch);
    };
    dfs(root);
}

static void emitContractedEdges(TreeNode* node,
                                const std::set<TreeNode*>& kept,
                                std::ofstream& out) {
    if (!node) return;
    for (auto* child : node->children) {
        TreeNode* cur = child;
        while (cur && kept.find(cur) == kept.end()) {
            if (cur->children.empty()) break;
            cur = cur->children[0];
        }
        if (cur) {
            const char *edge_color = (node->context.inversion == 1) ? "lightskyblue" : "black";
            out << "  node" << node->id
                << " -> node" << cur->id
                << " [color=\"" << edge_color << "\"];\n";
        }
    }
}

void writeCollapsedTreeDOT(TreeNode* root, const std::string& filename) {
    if (!root) return;

    std::set<TreeNode*> kept;
    collectKeptNodes(root, kept);

    std::ofstream out(filename.c_str());
    if (!out.is_open()) return;

    out << "digraph SiteTreeCollapsed {\n";
    out << "  node [shape=circle];\n";

    for (auto* n : kept) {
        out << "  node" << n->id
            << " [label=\""
            << n->id
            << "\\ninv=" << n->context.inversion
            << " pop=" << n->context.pop
            << "\\nt=" << n->time
            << "\"];\n";
    }

    for (auto* n : kept) {
        emitContractedEdges(n, kept, out);
    }

    out << "}\n";
}
