#ifndef SIMULATE_SMC_H
#define SIMULATE_SMC_H

#include <fstream>
#include <vector>
#include <string>
#include "treemod.h"
#include "parameters.h"

struct GeneFluxEvent_SMC {
    double endX = -1.0;
    unsigned long nodeId = 0;
};

struct SMCStepOutcome {
    std::vector<GeneFluxEvent_SMC> geneFluxEvents;
    std::vector<std::string> eventRows;
};

// Simulate unary migration events on the cut subtree until coalescence,
// then reattach to the main tree.
bool simulateSMCOnTree_SMC(TreeNode*& mainRoot,
                           TreeNode*& cutRoot,
                           double cutStartTime,
                           const Parameters::ParameterData& params,
                           const std::vector<std::vector<double>>& mig_prob,
                           const std::vector<std::pair<double, std::vector<std::vector<double>>>>& standard_mig_schedule,
                           const std::vector<std::pair<double, std::vector<std::vector<double>>>>& inverted_mig_schedule,
                           double currentHopX,
                           int hopIndex,
                           SMCStepOutcome* outcome,
                           int runIndex,
                           unsigned long& nextNodeId,
                           std::ofstream& evlog);

#endif
