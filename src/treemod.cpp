#include "treemod.h"
#include <unordered_map>
#include <fstream>
#include <limits>
#include "ran_mk.h"

static void buildRec(const SiteNode& src, TreeNode* parent, TreeNode* outNode,
                     std::unordered_map<unsigned long, TreeNode*>& idMap) {
    outNode->id = src.getNodeNumber();
    outNode->time = src.getTime();
    outNode->context = src.getContext();
    outNode->parent = parent;
    idMap[outNode->id] = outNode;
    for (const auto& childPtr : src.descendant) {
        if (!childPtr) continue;
        if (childPtr->getTime() < 0) continue;
        if (childPtr->getNodeNumber() == std::numeric_limits<unsigned long>::max()) continue;
        TreeNode* child = new TreeNode();
        outNode->children.push_back(child);
        buildRec(*childPtr, outNode, child, idMap);
    }
}

TreeNode* buildEditableTree(const SiteNode& root) {
    auto* out = new TreeNode();
    std::unordered_map<unsigned long, TreeNode*> idMap;
    buildRec(root, nullptr, out, idMap);
    return out;
}

static TreeNode* cloneRec(const TreeNode* src, TreeNode* parent) {
    if (!src) return nullptr;
    TreeNode* out = new TreeNode();
    out->id = src->id;
    out->time = src->time;
    out->context = src->context;
    out->parent = parent;
    for (auto* ch : src->children) {
        TreeNode* c = cloneRec(ch, out);
        if (c) out->children.push_back(c);
    }
    return out;
}

TreeNode* cloneTree(const TreeNode* root) {
    return cloneRec(root, nullptr);
}

static void freeRec(TreeNode* node) {
    if (!node) return;
    for (auto* ch : node->children) {
        freeRec(ch);
    }
    delete node;
}
void freeTree(TreeNode* root) {
    freeRec(root);
}

unsigned long getMaxId(TreeNode* root) {
    if (!root) return 0;
    unsigned long maxId = root->id;
    for (auto* ch : root->children) {
        unsigned long m = getMaxId(ch);
        if (m > maxId) maxId = m;
    }
    return maxId;
}

TreeNode* findNodeById(TreeNode* root, unsigned long id) {
    if (!root) return nullptr;
    if (root->id == id) return root;
    for (auto* ch : root->children) {
        if (auto* hit = findNodeById(ch, id)) return hit;
    }
    return nullptr;
}

TreeNode* addUnaryAbove(TreeNode* node, unsigned long newId, double newTime, const Context& newCtx) {
    if (!node) return nullptr;
    TreeNode* parent = node->parent;
    TreeNode* u = new TreeNode();
    u->id = newId;
    u->time = newTime;
    u->context = newCtx;
    u->parent = parent;
    u->children.push_back(node);

    if (parent) {
        for (auto& ch : parent->children) {
            if (ch == node) { ch = u; break; }
        }
    }
    node->parent = u;
    return u;
}

struct TimedEdgeCandidate {
    TreeNode* parent;
    TreeNode* child;
};

static void collectEdgesByTimeAndContext(TreeNode* node,
                                         double t,
                                         const Context& branchCtx,
                                         std::vector<TimedEdgeCandidate>& out) {
    if (!node) return;
    for (auto* ch : node->children) {
        if (node->time >= t && ch->time <= t && node->context == branchCtx) {
            TimedEdgeCandidate cand;
            cand.parent = node;
            cand.child = ch;
            out.push_back(cand);
        }
        collectEdgesByTimeAndContext(ch, t, branchCtx, out);
    }
}

bool reattachAtTimeWithContext(TreeNode*& mainRoot,
                               TreeNode* cutRoot,
                               double eventTime,
                               const Context& branchCtx,
                               unsigned long& nextId) {
    if (!mainRoot || !cutRoot) return false;

    std::vector<TimedEdgeCandidate> candidates;
    collectEdgesByTimeAndContext(mainRoot, eventTime, branchCtx, candidates);
    if (candidates.empty()) return false;

    const unsigned long pick = randint(0, static_cast<int>(candidates.size()) - 1);
    TreeNode* parent = candidates[pick].parent;
    TreeNode* child = candidates[pick].child;

    TreeNode* coal = new TreeNode();
    coal->id = nextId++;
    coal->time = eventTime;
    coal->context = branchCtx;
    coal->parent = parent;
    coal->children.push_back(child);
    coal->children.push_back(cutRoot);

    for (auto& ch : parent->children) {
        if (ch == child) {
            ch = coal;
            break;
        }
    }
    child->parent = coal;
    cutRoot->parent = coal;

    if (mainRoot == child && parent == nullptr) {
        mainRoot = coal;
    }
    return true;
}

static TreeNode* trimRetainedSiblingAfterCut(TreeNode* node, double cutTime) {
    while (node && node->children.size() == 1 && node->time > cutTime) {
        TreeNode* child = node->children.front();
        if (!(node->context == child->context)) break;
        node->children.clear();
        node->parent = nullptr;
        child->parent = nullptr;
        delete node;
        node = child;
    }
    return node;
}

bool cutAtNodeWithUnaryCleanup(TreeNode*& root, TreeNode* cutNode, double cutTime, TreeNode** outCutSubtree) {
    if (!root || !cutNode || !outCutSubtree) return false;
    if (!cutNode->parent) return false;

    // Track path from cut node upward through unary ancestors.
    TreeNode* pathTop = cutNode;
    TreeNode* curParent = cutNode->parent;
    std::vector<TreeNode*> unaryChain;
    while (curParent && curParent->children.size() == 1) {
        unaryChain.push_back(curParent);
        pathTop = curParent;
        curParent = curParent->parent;
    }

    TreeNode* coal = curParent; // first non-unary ancestor on this path
    TreeNode* coalParent = nullptr;
    TreeNode* sibling = nullptr;

    if (coal) {
        coalParent = coal->parent;
        for (auto* ch : coal->children) {
            if (ch != pathTop) {
                sibling = ch;
                break;
            }
        }
        sibling = trimRetainedSiblingAfterCut(sibling, cutTime);
    }

    // Detach cut subtree root.
    cutNode->parent = nullptr;
    *outCutSubtree = cutNode;

    // Rewire around coalescent node if present.
    if (coal && sibling) {
        if (coalParent) {
            for (auto& ch : coalParent->children) {
                if (ch == coal) {
                    ch = sibling;
                    break;
                }
            }
            sibling->parent = coalParent;
        } else {
            root = sibling;
            sibling->parent = nullptr;
        }
    }
    if (!coal && !unaryChain.empty() && unaryChain.back() == root) {
        root = nullptr;
    }

    // Remove and delete unary chain nodes between cutNode and coalescent.
    for (auto* u : unaryChain) {
        u->children.clear();
        u->parent = nullptr;
        delete u;
    }

    // Delete coalescent node itself as requested after rewiring sibling.
    if (coal) {
        coal->children.clear();
        coal->parent = nullptr;
        delete coal;
    }

    return true;
}

bool cutEdgeRandomWithCleanup(TreeNode*& root,
                              TreeNode* parent,
                              TreeNode* child,
                              TreeNode** outCutSubtree,
                              double* outCutTime) {
    if (!root || !parent || !child || !outCutSubtree) return false;
    if (child->parent != parent) return false;
    const double t_parent = parent->time;
    const double t_child = child->time;
    if (t_parent <= t_child) return false;

    const double t_cut = randreal(t_child, t_parent);
    if (outCutTime) {
        *outCutTime = t_cut;
    }

    return cutAtNodeWithUnaryCleanup(root, child, t_cut, outCutSubtree);
}

void trimUnaryRootStem(TreeNode*& root) {
    while (root && root->children.size() == 1) {
        TreeNode* oldRoot = root;
        TreeNode* newRoot = root->children.front();
        oldRoot->children.clear();
        oldRoot->parent = nullptr;
        newRoot->parent = nullptr;
        root = newRoot;
        delete oldRoot;
    }
}

static void gatherNodes(TreeNode* root, std::vector<TreeNode*>& out) {
    if (!root) return;
    out.push_back(root);
    for (auto* ch : root->children) gatherNodes(ch, out);
}

void writeTreeDOT(TreeNode* root, const std::string& filename) {
    std::vector<TreeNode*> nodes;
    gatherNodes(root, nodes);

    std::ofstream outFile(filename.c_str());
    if (!outFile.is_open()) return;

    outFile << "digraph SiteTreeModified {\n";
    outFile << "  node [shape=circle];\n";

    for (auto* node : nodes) {
        outFile << "  node" << node->id
                << " [label=\""
                << node->id
                << "\\ninv=" << node->context.inversion
                << " pop=" << node->context.pop
                << "\\nt=" << node->time
                << "\"];\n";
    }

    for (auto* node : nodes) {
        for (auto* ch : node->children) {
            const char *edge_color = (node->context.inversion == 1) ? "lightskyblue" : "black";
            outFile << "  node" << node->id
                    << " -> node" << ch->id
                    << " [color=\"" << edge_color << "\"];\n";
        }
    }

    outFile << "}\n";
}
