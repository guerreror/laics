#ifndef SMC_HELPERS_H
#define SMC_HELPERS_H

#include <string>
#include "treemod.h"

// Writes a visualization-only collapsed tree:
// unary chains are contracted, while coalescent/branching nodes and leaves remain.
// This does not mutate the input tree.
void writeCollapsedTreeDOT(TreeNode* root, const std::string& filename);

#endif

