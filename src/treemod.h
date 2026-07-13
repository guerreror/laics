#ifndef TREEMOD_H
#define TREEMOD_H

#include <vector>
#include "typedefs.h"
#include "sitenode.h"

struct TreeNode {
    unsigned long id;
    double time;
    Context context;
    TreeNode* parent = nullptr;
    std::vector<TreeNode*> children;
};

TreeNode* buildEditableTree(const SiteNode& root);
TreeNode* cloneTree(const TreeNode* root);
void freeTree(TreeNode* root);

TreeNode* findNodeById(TreeNode* root, unsigned long id);
bool cutAtNodeWithUnaryCleanup(TreeNode*& root, TreeNode* cutNode, TreeNode** outCutSubtree);
bool cutEdgeRandomWithCleanup(TreeNode*& root,
                              TreeNode* parent,
                              TreeNode* child,
                              TreeNode** outCutSubtree,
                              double* outCutTime);
void trimUnaryRootStem(TreeNode*& root);

unsigned long getMaxId(TreeNode* root);
TreeNode* addUnaryAbove(TreeNode* node, unsigned long newId, double newTime, const Context& newCtx);
bool reattachAtTimeWithContext(TreeNode*& mainRoot,
                               TreeNode* cutRoot,
                               double eventTime,
                               const Context& branchCtx,
                               unsigned long& nextId);

void writeTreeDOT(TreeNode* root, const std::string& filename);

#endif
