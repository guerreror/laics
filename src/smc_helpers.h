#ifndef SMC_HELPERS_H
#define SMC_HELPERS_H

#include <string>
#include <memory>
#include "argnode.h"
#include "treemod.h"

// Writes a visualization-only collapsed tree:
// unary chains are contracted, while coalescent/branching nodes and leaves remain.
// This does not mutate the input tree.
void writeCollapsedTreeDOT(TreeNode* root, const std::string& filename);

// Builds the local x0 tree from an ARG while preserving unary ARG event nodes
// on branches that carry x0. The returned tree is owned by the caller.
TreeNode* buildX0TreeFromARGPreserveUnary(double x0, std::shared_ptr<ARGNode> argRoot);

#endif
