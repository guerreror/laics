#ifndef SIMULATE_SMC_H
#define SIMULATE_SMC_H

#include <vector>
#include <string>
#include "treemod.h"
#include "parameters.h"

struct GeneFluxEvent_SMC {
    double startX = -1.0;
    double endX = -1.0;
    unsigned long nodeId = 0;
    std::string type;
};

struct SMCEpochs_SMC {
    std::vector<double> breaks;
};

struct SMCStepOutcome {
    bool coalesced = false;
    bool hitRootLimit = false;
    double stopTime = 0.0;
    std::vector<GeneFluxEvent_SMC> geneFluxEvents;
    std::vector<std::string> eventRows;
};

SMCEpochs_SMC buildEpochBreaks_SMC(TreeNode* mainTree, double cut_time);

// Simulate unary migration events on the cut subtree until coalescence,
// then reattach to the main tree.
bool simulateSMCOnTree_SMC(TreeNode*& mainRoot,
                           TreeNode*& cutRoot,
                           double cutStartTime,
                           const Parameters::ParameterData& params,
                           const std::vector<std::vector<double>>& mig_prob,
                           double currentHopX,
                           int hopIndex,
                           SMCStepOutcome* outcome,
                           const std::string& eventLogPath);

#endif
