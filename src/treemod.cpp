#include "treemod.h"
#include "ran_mk.h"

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

static void noteDeleted(TreeNode* node, std::vector<unsigned long>* deletedIds) {
    if (node && deletedIds) deletedIds->push_back(node->id);
}

static TreeNode* trimRetainedSiblingAfterCut(TreeNode* node,
                                             double cutTime,
                                             std::vector<unsigned long>* deletedIds) {
    while (node && node->children.size() == 1 && node->time > cutTime) {
        TreeNode* child = node->children.front();
        if (!(node->context == child->context)) break;
        node->children.clear();
        node->parent = nullptr;
        child->parent = nullptr;
        noteDeleted(node, deletedIds);
        delete node;
        node = child;
    }
    return node;
}

static bool cutAtNodeWithUnaryCleanupImpl(TreeNode*& root,
                                          TreeNode* cutNode,
                                          double cutTime,
                                          TreeNode** outCutSubtree,
                                          std::vector<unsigned long>* deletedIds) {
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
        sibling = trimRetainedSiblingAfterCut(sibling, cutTime, deletedIds);
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
        noteDeleted(u, deletedIds);
        delete u;
    }

    // Delete coalescent node itself as requested after rewiring sibling.
    if (coal) {
        coal->children.clear();
        coal->parent = nullptr;
        noteDeleted(coal, deletedIds);
        delete coal;
    }

    return true;
}

bool cutAtNodeWithUnaryCleanupCollect(TreeNode*& root,
                                      TreeNode* cutNode,
                                      double cutTime,
                                      TreeNode** outCutSubtree,
                                      std::vector<unsigned long>& deletedIds) {
    return cutAtNodeWithUnaryCleanupImpl(root, cutNode, cutTime, outCutSubtree, &deletedIds);
}

bool cutEdgeRandomWithCleanupCollect(TreeNode*& root,
                                     TreeNode* parent,
                                     TreeNode* child,
                                     TreeNode** outCutSubtree,
                                     double* outCutTime,
                                     std::vector<unsigned long>& deletedIds) {
    if (!root || !parent || !child || !outCutSubtree) return false;
    if (child->parent != parent) return false;
    const double t_parent = parent->time;
    const double t_child = child->time;
    if (t_parent <= t_child) return false;

    const double t_cut = randreal(t_child, t_parent);
    if (outCutTime) {
        *outCutTime = t_cut;
    }

    return cutAtNodeWithUnaryCleanupCollect(root, child, t_cut, outCutSubtree, deletedIds);
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
