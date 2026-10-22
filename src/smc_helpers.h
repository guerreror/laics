#ifndef SMC_HELPERS_H
#define SMC_HELPERS_H

#include <string>
#include <memory>
#include <vector>
#include "argnode.h"
#include "parameters.h"
#include "treemod.h"

struct SMCActiveState {
    std::vector<unsigned int> popSizes;
    std::vector<double> invFreqs;
};

const SMCActiveState& activeStateAtTime_SMC(const Parameters::ParameterData& params, double t);

// Builds the local x0 tree from an ARG while preserving unary ARG event nodes
// on branches that carry x0. The returned tree is owned by the caller.
TreeNode* buildX0TreeFromARGPreserveUnary(double x0, std::shared_ptr<ARGNode> argRoot);

#endif
