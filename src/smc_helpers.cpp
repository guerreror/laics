#include "smc_helpers.h"

#include <fstream>
#include <set>
#include <vector>
#include <functional>

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

