#ifndef TREEMOD_H
#define TREEMOD_H

#include <vector>
#include "typedefs.h"

struct TreeNode {
    unsigned long id;
    double time;
    Context context;
    TreeNode* parent = nullptr;
    std::vector<TreeNode*> children;
};

void freeTree(TreeNode* root);

TreeNode* findNodeById(TreeNode* root, unsigned long id);
bool cutAtNodeWithUnaryCleanupCollect(TreeNode*& root,
                                      TreeNode* cutNode,
                                      double cutTime,
                                      TreeNode** outCutSubtree,
                                      std::vector<unsigned long>& deletedIds);
bool cutEdgeRandomWithCleanupCollect(TreeNode*& root,
                                     TreeNode* parent,
                                     TreeNode* child,
                                     TreeNode** outCutSubtree,
                                     double* outCutTime,
                                     std::vector<unsigned long>& deletedIds);
void trimUnaryRootStem(TreeNode*& root);

unsigned long getMaxId(TreeNode* root);
TreeNode* addUnaryAbove(TreeNode* node, unsigned long newId, double newTime, const Context& newCtx);

#endif
