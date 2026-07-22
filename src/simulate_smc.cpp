/** simulate_smc.cpp
 *  SMC-specific simulation: migration + coalescence + gene flux on editable trees.
 */

#include "simulate_smc.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <utility>
#include <vector>

#include <boost/math/special_functions/binomial.hpp>

#include "argnode.h"
#include "chromosome.h"
#include "ran_mk.h"
#include "world.h"

using boost::math::binomial_coefficient;
using std::map;
using std::vector;

static const double SMC_DEBUG_MIGRATION_BOOST = 1.0;

static void recordEventRow_SMC(SMCStepOutcome* outcome,
                               std::ofstream& evlog,
                               int hopIndex,
                               double currentHopX,
                               const std::string& eventName,
                               double eventTime) {
    std::ostringstream row;
    row << hopIndex << "," << currentHopX << "," << eventName << ","
        << eventTime << ",,,\n";
    if (outcome) {
        outcome->eventRows.push_back(row.str());
    }
    if (evlog.is_open()) {
        evlog << row.str();
    }
}

static void gatherTimes_SMC(TreeNode* node, std::vector<double>& out) {
    if (!node) return;
    out.push_back(node->time);
    for (auto* ch : node->children) gatherTimes_SMC(ch, out);
}

SMCEpochs_SMC buildEpochBreaks_SMC(TreeNode* mainTree, double cut_time) {
    std::vector<double> times;
    gatherTimes_SMC(mainTree, times);
    SMCEpochs_SMC epochs;
    for (double t : times) {
        if (t > cut_time) epochs.breaks.push_back(t);
    }
    std::sort(epochs.breaks.begin(), epochs.breaks.end());
    epochs.breaks.erase(std::unique(epochs.breaks.begin(), epochs.breaks.end()),
                        epochs.breaks.end());
    return epochs;
}

static double computeTotalM_SMC(TreeNode* cutRoot,
                                const Parameters::ParameterData& params,
                                const std::vector<std::vector<double>>& mig_prob) {
    if (!cutRoot) return 0.0;
    const unsigned int pop = cutRoot->context.pop;
    const unsigned int nPops = params.popSizeVec.size();
    if (pop >= nPops) return 0.0;

    double b_mig_total = 0.0;
    for (unsigned int l = 0; l < nPops; ++l) {
        if (l != pop) b_mig_total += mig_prob.at(pop).at(l);
    }
    return b_mig_total;
}

static double computeTotalC_SMC(const Parameters::ParameterData& params) {
    if (params.totalPopSize == 0) return 0.0;
    return 1.0 / static_cast<double>(params.totalPopSize);
}

static double computeTotalG_SMC(const TreeNode* cutRoot,
                                const Parameters::ParameterData& params,
                                double /*currentHopX*/) {
    if (!cutRoot) return 0.0;
    const unsigned int pop = cutRoot->context.pop;
    if (pop >= params.initialFreqs.size()) return 0.0;

    const double invFreq = params.initialFreqs[pop];
    const double stdFreq = 1.0 - invFreq;
    if (cutRoot->context.inversion == 0) {
        return params.phi * invFreq;
    }
    return params.phi * stdFreq;
}

static double drawGeneFluxSegmentLength_SMC(double /*currentHopX*/) {
    return 0.01;
}

static double triangleHeightAtX_SMC(double x, const Segment& range, double peakHeight) {
    const double L = range.L;
    const double R = range.R;
    if (R <= L || x <= L || x >= R) return 0.0;

    const double mid = 0.5 * (L + R);
    if (x <= mid) {
        return peakHeight * ((x - L) / (mid - L));
    }
    return peakHeight * ((R - x) / (R - mid));
}

static std::string pickGeneFluxType_SMC(double x, const Parameters::ParameterData& params) {
    const double gcHeight = std::max(0.0, params.gcRate);
    const double drHeight = triangleHeightAtX_SMC(x, params.smcRange, std::max(0.0, params.drRate));
    const double totalHeight = gcHeight + drHeight;
    if (totalHeight <= 0.0) return "GC";

    const double u = randreal(0.0, totalHeight);
    return (u < gcHeight) ? "GC" : "DR";
}

static double drawDoubleRecombinationEnd_SMC(double startX, const Segment& range) {
    const double mid = 0.5 * (range.L + range.R);
    const double lower = std::max(startX, mid);
    const double upper = range.R;
    if (upper <= lower) return upper;
    return randreal(lower, upper);
}

static GeneFluxEvent_SMC makeGeneFluxSegment_SMC(double startX,
                                                 unsigned long nodeId,
                                                 const std::string& type,
                                                 const Parameters::ParameterData& params) {
    GeneFluxEvent_SMC evt;
    evt.startX = std::max(params.smcRange.L, std::min(params.smcRange.R, startX));
    evt.nodeId = nodeId;
    evt.type = type;

    if (type == "DR") {
        evt.endX = drawDoubleRecombinationEnd_SMC(evt.startX, params.smcRange);
    } else {
        const double J = drawGeneFluxSegmentLength_SMC(evt.startX);
        evt.endX = evt.startX + J;
    }
    evt.endX = std::max(evt.startX, std::min(params.smcRange.R, evt.endX));
    return evt;
}

static unsigned int pickMigrationDest_SMC(unsigned int fromPop,
                                          const std::vector<std::vector<double>>& mig_prob) {
    const auto& row = mig_prob.at(fromPop);
    double roll = randreal(0, 1);
    double cum = 0.0;
    for (unsigned int i = 0; i < row.size(); ++i) {
        cum += row[i];
        if (roll <= cum) return i;
    }
    return fromPop;
}

static bool canCoalesceByContext_SMC(const TreeNode* a, const TreeNode* b) {
    if (!a || !b) return false;
    return a->context == b->context;
}

static Context contextAfterSpeciationEvents_SMC(Context ctx,
                                                const Parameters::ParameterData& params,
                                                double t) {
    if (params.speciation.empty() || params.speciation[0] != 1) {
        return ctx;
    }

    const double eps = 1e-9;
    for (size_t i = 1; i + 4 < params.speciation.size(); i += 5) {
        const unsigned int sink = static_cast<unsigned int>(params.speciation[i]);
        const unsigned int source = static_cast<unsigned int>(params.speciation[i + 1]);
        const double eventTime = params.speciation[i + 2];
        if (t + eps < eventTime) {
            continue;
        }

        const unsigned int newSink = (sink > source) ? (sink - 1) : sink;
        if (ctx.pop == source) {
            ctx.pop = newSink;
        } else if (ctx.pop > source) {
            ctx.pop -= 1;
        }
    }
    return ctx;
}

static Context contextAfterEpochEvents_SMC(Context ctx,
                                           const Parameters::ParameterData& params,
                                           double t) {
    ctx = contextAfterSpeciationEvents_SMC(ctx, params, t);

    const double eps = 1e-9;
    if (params.inv_age > 0 && t + eps >= static_cast<double>(params.inv_age) &&
        ctx.inversion == 1) {
        // This mirrors World::freqStepToLoss(): beyond inversion age all
        // inversion contexts collapse to the ancestral origin context.
        ctx.pop = 0;
        ctx.inversion = 0;
    }
    return ctx;
}

static bool nextModelEpochAfter_SMC(const Parameters::ParameterData& params,
                                    double currentTime,
                                    double& nextEpoch) {
    const double eps = 1e-9;
    bool found = false;
    double best = 0.0;
    auto consider = [&](double t) {
        if (t > currentTime + eps && (!found || t < best)) {
            best = t;
            found = true;
        }
    };

    if (params.inv_age > 0) {
        consider(static_cast<double>(params.inv_age));
    }
    if (!params.speciation.empty() && params.speciation[0] == 1) {
        for (size_t i = 1; i + 4 < params.speciation.size(); i += 5) {
            consider(params.speciation[i + 2]);
        }
    }

    if (found) {
        nextEpoch = best;
    }
    return found;
}

static void applyEpochEventsToLineageAtTime_SMC(TreeNode*& lineage,
                                                double t,
                                                const Parameters::ParameterData& params,
                                                unsigned long& nextId) {
    if (!lineage) return;
    const Context newCtx = contextAfterEpochEvents_SMC(lineage->context, params, t);
    if (newCtx == lineage->context) {
        return;
    }

    TreeNode* nr = addUnaryAbove(lineage, nextId++, t, newCtx);
    if (nr && nr->parent == nullptr) {
        lineage = nr;
    }
}

static bool resolveAboveRootByMiniSMC_SMC(TreeNode*& mainRoot,
                                          TreeNode*& cutRoot,
                                          double cutLineageTime,
                                          const Parameters::ParameterData& params,
                                          const std::vector<std::vector<double>>& mig_prob,
                                          double currentHopX,
                                          SMCStepOutcome* outcome,
                                          std::ofstream& evlog,
                                          int hopIndex,
                                          double root_time) {
    if (!mainRoot || !cutRoot) {
        std::cerr << "SMC fallback failed: missing main or cut lineage.\n";
        return false;
    }

    unsigned long nextId = std::max(getMaxId(mainRoot), getMaxId(cutRoot)) + 1;
    double currentTime = std::max(root_time, cutLineageTime);
    if (mainRoot->time < currentTime) {
        TreeNode* nr = addUnaryAbove(mainRoot, nextId++, currentTime, mainRoot->context);
        if (nr && nr->parent == nullptr) mainRoot = nr;
    }
    if (cutRoot->time < currentTime) {
        TreeNode* nr = addUnaryAbove(cutRoot, nextId++, currentTime, cutRoot->context);
        if (nr && nr->parent == nullptr) cutRoot = nr;
    }
    applyEpochEventsToLineageAtTime_SMC(mainRoot, currentTime, params, nextId);
    applyEpochEventsToLineageAtTime_SMC(cutRoot, currentTime, params, nextId);

    std::vector<TreeNode*> lineages;
    lineages.push_back(mainRoot);
    lineages.push_back(cutRoot);

    while (lineages.size() > 1) {
        std::vector<double> totalM_i(lineages.size(), 0.0);
        std::vector<double> totalG_i(lineages.size(), 0.0);
        double totalM = 0.0;
        double totalG = 0.0;
        for (size_t i = 0; i < lineages.size(); ++i) {
            totalM_i[i] = computeTotalM_SMC(lineages[i], params, mig_prob) * SMC_DEBUG_MIGRATION_BOOST;
            totalG_i[i] = computeTotalG_SMC(lineages[i], params, currentHopX);
            totalM += totalM_i[i];
            totalG += totalG_i[i];
        }
        const double totalC = computeTotalC_SMC(params);
        const double Rate = totalC + totalM + totalG;
        if (Rate <= 0.0) {
            std::cerr << "SMC fallback failed: total event rate is zero.\n";
            return false;
        }

        const double startTime = currentTime;
        const double dt = randexp(Rate);
        const double event_time = currentTime + dt;
        double nextEpoch = 0.0;
        if (nextModelEpochAfter_SMC(params, startTime, nextEpoch) &&
            event_time >= nextEpoch) {
            currentTime = nextEpoch;
            for (auto*& lineage : lineages) {
                applyEpochEventsToLineageAtTime_SMC(lineage, currentTime, params, nextId);
            }
            continue;
        }

        currentTime = event_time;
        for (auto*& lineage : lineages) {
            applyEpochEventsToLineageAtTime_SMC(lineage, currentTime, params, nextId);
        }

        const double roll = randreal(0, Rate);
        if (roll < totalC) {
            std::vector<std::pair<size_t, size_t>> pairs;
            for (size_t i = 0; i < lineages.size(); ++i) {
                for (size_t j = i + 1; j < lineages.size(); ++j) {
                    if (canCoalesceByContext_SMC(lineages[i], lineages[j])) {
                        pairs.push_back({i, j});
                    }
                }
            }
            if (pairs.empty()) {
                if (totalM <= 0.0 && totalG <= 0.0) {
                    double nextEpoch = 0.0;
                    if (nextModelEpochAfter_SMC(params, currentTime, nextEpoch)) {
                        currentTime = nextEpoch;
                        for (auto*& lineage : lineages) {
                            applyEpochEventsToLineageAtTime_SMC(lineage, currentTime, params, nextId);
                        }
                        continue;
                    }
                    std::ostringstream reason;
                    reason << "incompatible_context_fallback"
                           << "_main_p" << (lineages.size() > 0 ? lineages[0]->context.pop : 9999)
                           << "_i" << (lineages.size() > 0 ? lineages[0]->context.inversion : 9999)
                           << "_cut_p" << (lineages.size() > 1 ? lineages[1]->context.pop : 9999)
                           << "_i" << (lineages.size() > 1 ? lineages[1]->context.inversion : 9999);
                    recordEventRow_SMC(outcome, evlog, hopIndex, currentHopX,
                                       reason.str(), currentTime);
                    return false;
                }
                continue;
            }

            const auto pick = pairs[randint(0, static_cast<int>(pairs.size()) - 1)];
            const size_t i = pick.first;
            const size_t j = pick.second;
            TreeNode* coal = new TreeNode();
            coal->id = nextId++;
            coal->time = event_time;
            coal->context = lineages[i]->context;
            coal->children.push_back(lineages[i]);
            coal->children.push_back(lineages[j]);
            lineages[i]->parent = coal;
            lineages[j]->parent = coal;

            lineages.erase(lineages.begin() + static_cast<long>(j));
            lineages.erase(lineages.begin() + static_cast<long>(i));
            lineages.push_back(coal);
            recordEventRow_SMC(outcome, evlog, hopIndex, currentHopX, "coalescence_fallback", event_time);
            continue;
        }

        double rem = roll - totalC;
        bool handled = false;
        for (size_t i = 0; i < lineages.size(); ++i) {
            if (rem < totalM_i[i]) {
                Context newCtx = lineages[i]->context;
                newCtx.pop = pickMigrationDest_SMC(newCtx.pop, mig_prob);
                TreeNode* nr = addUnaryAbove(lineages[i], nextId++, event_time, newCtx);
                if (nr && nr->parent == nullptr) lineages[i] = nr;
                recordEventRow_SMC(outcome, evlog, hopIndex, currentHopX, "migration_fallback", event_time);
                handled = true;
                break;
            }
            rem -= totalM_i[i];

            if (rem < totalG_i[i]) {
                Context newCtx = lineages[i]->context;
                newCtx.inversion = (newCtx.inversion == 0 ? 1 : 0);
                TreeNode* nr = addUnaryAbove(lineages[i], nextId++, event_time, newCtx);
                if (nr && nr->parent == nullptr) lineages[i] = nr;
                const std::string geneFluxType = pickGeneFluxType_SMC(currentHopX, params);
                if (outcome) {
                    GeneFluxEvent_SMC evt = makeGeneFluxSegment_SMC(currentHopX,
                                                                    nr ? nr->id : 0,
                                                                    geneFluxType,
                                                                    params);
                    outcome->geneFluxEvents.push_back(evt);
                }
                recordEventRow_SMC(outcome, evlog, hopIndex, currentHopX,
                                   "gene_flux_fallback_" + geneFluxType, event_time);
                handled = true;
                break;
            }
            rem -= totalG_i[i];
        }
        if (!handled) {
            std::cerr << "SMC fallback failed: sampled event was outside all rate channels.\n";
            return false;
        }
    }

    mainRoot = lineages[0];
    if (outcome) {
        outcome->coalesced = true;
        outcome->hitRootLimit = false;
        outcome->stopTime = mainRoot->time;
    }
    return true;
}

bool simulateSMCOnTree_SMC(TreeNode*& mainRoot,
                           TreeNode*& cutRoot,
                           double cutStartTime,
                           const Parameters::ParameterData& params,
                           const std::vector<std::vector<double>>& mig_prob,
                           double currentHopX,
                           int hopIndex,
                           SMCStepOutcome* outcome,
                           const std::string& eventLogPath) {
    if (!mainRoot || !cutRoot) {
        std::cerr << "SMC reattachment failed: missing main or cut lineage.\n";
        return false;
    }
    if (outcome) {
        *outcome = SMCStepOutcome{};
    }

    const double root_time = mainRoot->time;
    std::ofstream evlog;
    if (!eventLogPath.empty()) {
        evlog.open(eventLogPath.c_str(), std::ios::app);
    }
    if (evlog.is_open() && evlog.tellp() == 0) {
        evlog << "hop,current_x,event,event_time,raw_delta_x,used_delta_x,next_x\n";
    }
    double nextGeneFluxStartX = currentHopX;
    double lineageTime = std::max(cutStartTime, cutRoot->time);

    while (true) {
        SMCEpochs_SMC epochs = buildEpochBreaks_SMC(mainRoot, lineageTime);

        double totalM = computeTotalM_SMC(cutRoot, params, mig_prob) * SMC_DEBUG_MIGRATION_BOOST;
        double totalC = computeTotalC_SMC(params);
        double totalG = computeTotalG_SMC(cutRoot, params, currentHopX);
        double Rate = totalM + totalC + totalG;
        if (Rate <= 0.0) {
            std::cerr << "SMC reattachment failed: total event rate is zero.\n";
            return false;
        }

        double waiting_t = randexp(Rate);
        double event_time = lineageTime + waiting_t;

        bool crossedEpoch = false;
        double nextEpoch = 0.0;
        if (!epochs.breaks.empty()) {
            auto it = std::upper_bound(epochs.breaks.begin(), epochs.breaks.end(), lineageTime);
            if (it != epochs.breaks.end() && event_time >= *it) {
                nextEpoch = *it;
                crossedEpoch = true;
            }
        }
        double modelEpoch = 0.0;
        if (nextModelEpochAfter_SMC(params, lineageTime, modelEpoch) &&
            event_time >= modelEpoch &&
            (!crossedEpoch || modelEpoch < nextEpoch)) {
            nextEpoch = modelEpoch;
            crossedEpoch = true;
        }
        if (crossedEpoch && nextEpoch < root_time) {
            lineageTime = nextEpoch;
            unsigned long nextId = std::max(getMaxId(mainRoot), getMaxId(cutRoot)) + 1;
            applyEpochEventsToLineageAtTime_SMC(cutRoot, lineageTime, params, nextId);
            continue;
        }

        if (event_time >= root_time) {
            if (resolveAboveRootByMiniSMC_SMC(mainRoot, cutRoot, lineageTime, params, mig_prob, currentHopX, outcome, evlog, hopIndex, root_time)) {
                return true;
            }
            std::cerr << "SMC reattachment failed: above-root simulation did not coalesce.\n";
            if (outcome) {
                outcome->coalesced = false;
                outcome->hitRootLimit = true;
                outcome->stopTime = root_time;
            }
            return false;
        }

        double roll = randreal(0, Rate);
        bool is_migration = (roll < totalM);
        bool is_coalescence = (!is_migration && roll < (totalM + totalC));
        bool is_gene_flux = (!is_migration && !is_coalescence);

        if (is_migration) {
            unsigned long nextId = std::max(getMaxId(mainRoot), getMaxId(cutRoot)) + 1;
            Context newCtx = cutRoot->context;
            unsigned int dest = pickMigrationDest_SMC(newCtx.pop, mig_prob);
            newCtx.pop = dest;
            TreeNode* newRoot = addUnaryAbove(cutRoot, nextId, event_time, newCtx);
            if (newRoot->parent == nullptr) {
                cutRoot = newRoot;
            }
            lineageTime = event_time;
            recordEventRow_SMC(outcome, evlog, hopIndex, currentHopX, "migration", event_time);
        } else if (is_coalescence) {
            unsigned long nextId = std::max(getMaxId(mainRoot), getMaxId(cutRoot)) + 1;
            if (reattachAtTimeWithContext(mainRoot, cutRoot, event_time, cutRoot->context, nextId)) {
                recordEventRow_SMC(outcome, evlog, hopIndex, currentHopX, "coalescence", event_time);
                if (outcome) {
                    outcome->coalesced = true;
                    outcome->hitRootLimit = false;
                    outcome->stopTime = event_time;
                }
                return true;
            }

            if (!epochs.breaks.empty()) {
                auto it = std::upper_bound(epochs.breaks.begin(), epochs.breaks.end(), lineageTime);
                if (it != epochs.breaks.end()) {
                    lineageTime = *it;
                    unsigned long nextId = std::max(getMaxId(mainRoot), getMaxId(cutRoot)) + 1;
                    applyEpochEventsToLineageAtTime_SMC(cutRoot, lineageTime, params, nextId);
                    continue;
                }
            }
            if (outcome) {
                outcome->coalesced = false;
                outcome->hitRootLimit = false;
                outcome->stopTime = lineageTime;
            }
            std::cerr << "SMC reattachment failed: no compatible main-tree edge and no later epoch.\n";
            return false;
        } else if (is_gene_flux) {
            unsigned long nextId = std::max(getMaxId(mainRoot), getMaxId(cutRoot)) + 1;
            Context newCtx = cutRoot->context;
            newCtx.inversion = (newCtx.inversion == 0 ? 1 : 0);
            TreeNode* newRoot = addUnaryAbove(cutRoot, nextId, event_time, newCtx);
            if (newRoot->parent == nullptr) {
                cutRoot = newRoot;
            }
            lineageTime = event_time;

            const std::string geneFluxType = pickGeneFluxType_SMC(nextGeneFluxStartX, params);
            GeneFluxEvent_SMC evt = makeGeneFluxSegment_SMC(nextGeneFluxStartX,
                                                            nextId,
                                                            geneFluxType,
                                                            params);
            if (outcome) {
                outcome->geneFluxEvents.push_back(evt);
            }
            nextGeneFluxStartX = evt.endX;
            recordEventRow_SMC(outcome, evlog, hopIndex, currentHopX,
                               "gene_flux_" + geneFluxType, event_time);
        }
    }

}

unsigned short World::simulateGeneration_SMC(vector<vector<double>>& mig_prob) {
    DBG(">> Simulating Generation (SMC) " << worldData->generation);

    vector<vector<double>> migmap;
    vector<double> mRate;
    vector<double> cRate;
    double totalM = 0;
    double totalC = 0;

    for (cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i) {
        Context cxt = i->first;
        unsigned int k = ctxtNcarriers(cxt);
        if (k > 0) {
            unsigned long clustID = i->second;
            int pop = i->first.pop;
            double clustSize = worldData->freq.at(clustID) * worldData->popSize.at(pop);

            vector<double> mfreq = getSizebyPop(cxt);
            vector<double> migs;
            double b_mig_total = 0;
            for (int l = 0; l < worldData->nPops; ++l) {
                if (l != pop) {
                    double q1q2 = 0;
                    if (mfreq.at(pop) > 0) q1q2 = mfreq.at(l) / mfreq.at(pop);
                    b_mig_total += mig_prob.at(pop).at(l) * q1q2;
                    migs.push_back(b_mig_total);
                } else {
                    migs.push_back(0);
                }
            }
            migmap.push_back(migs);

            mRate.push_back(k * b_mig_total);
            totalM += mRate[clustID];

            if (k > 1) {
                unsigned int a = 2;
                if (clustSize == 0) {
                    std::cerr << "Error in simulateGeneration_SMC(): Carriers in empty context\n";
                }
                double ka = binomial_coefficient<double>(k, a) / clustSize;
                cRate.push_back(ka);
                totalC += ka;
            }
        } else {
            vector<double> dum(1, 0);
            migmap.push_back(dum);
            mRate.push_back(0);
        }
    }

    int nEvents = 0;
    double Rate = totalM + totalC;
    if (Rate == 0.0) {
        if (!worldData->epochs_over) {
            updateToNextEpoch_SMC();
            return 0;
        }
        forceAllCoal();
        return 0;
    }
    double waiting_t = randexp(Rate);

    if (timeCheck_SMC(waiting_t, Rate)) {
        updateToNextEpoch_SMC();
    } else {
        worldData->generation += waiting_t;
        double event = randreal(0, 1);
        if (Rate > 0) {
            if (event < totalM / Rate) nEvents += migrateEvent_SMC(migmap, mRate, totalM);
            else nEvents += coalesceEvent_SMC(cRate, totalC);
        }
    }

    if (sitesCoalesced()) {
        if (!worldData->epochs_over) {
            updateToNextEpoch_SMC();
        } else {
            forceAllCoal();
        }
    }

    return nEvents;
}

bool World::timeCheck_SMC(double waiting, double Rate) {
    double t = waiting + worldData->generation;
    if (!worldData->epochs_over) {
        if (t >= worldData->epoch_breaks.at(worldData->current_epoch)) {
            return true;
        }
        if (Rate == 0) {
            return true;
        }
    } else if (Rate == 0) {
        return false;
    }
    return false;
}

void World::updateToNextEpoch_SMC() {
    int outof_epoch = worldData->current_epoch;
    int into_epoch = worldData->current_epoch + 1;
    worldData->generation = worldData->epoch_breaks[outof_epoch];

    if (into_epoch <= static_cast<int>(worldData->epoch_breaks.size())) {
        switch (worldData->epochType[outof_epoch]) {
            case 0:
                freqStepToLoss();
                break;
            case 1:
                speciation();
                break;
            case 2:
                demoChange();
                break;
            default:
                break;
        }
    } else {
        std::cerr << "Error in World::updateToNextEpoch_SMC(). Current epoch > epoch vector size\n";
        exit(1);
    }

    if (into_epoch == static_cast<int>(worldData->epoch_breaks.size())) {
        worldData->epochs_over = true;
    }
    worldData->current_epoch = into_epoch;
}

unsigned short World::migrateEvent_SMC(vector<vector<double>>& mig_prob, vector<double>& rate, double total) {
    map<double, int> whichClust;
    double add = 0;
    for (cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i) {
        if (ctxtNcarriers(i->first) != 0) {
            whichClust[(rate[i->second] + add) / total] = i->second;
            add += rate[i->second];
        }
    }

    unsigned long c = whichClust.lower_bound(randreal(0, 1))->second;
    unsigned long who = randint(0, worldData->carriers->at(c).size() - 1);
    std::shared_ptr<Chromosome> chrom = worldData->carriers->at(c).at(who);

    vector<std::shared_ptr<Chromosome>>::iterator pos = worldData->carriers->at(c).begin() + who;
    worldData->carriers->at(c).erase(pos);

    int whereto = 0;
    if (chrom->getContext().pop == 0) whereto = 1; // current SMC path assumes two extant pops
    chrom->setPopulation(whereto);

    int newC = cluster[chrom->getContext()];
    worldData->carriers->at(newC).push_back(chrom);

    std::shared_ptr<ARGNode> newNode;
    newNode.reset(new ARGNode(worldData->argNodeVec.size(), chrom, worldData->generation));
    worldData->argNodeVec.push_back(newNode);
    chrom->setDescendant(newNode);

    return 1;
}

unsigned short World::coalesceEvent_SMC(vector<double>& rate, double total) {
    map<double, int> whichClust;

    int id = 0;
    double add = 0;
    for (cluster_t::iterator i = cluster.begin(); i != cluster.end(); ++i) {
        if (ctxtNcarriers(i->first) > 1) {
            whichClust[(rate.at(id) + add) / total] = i->second;
            add += rate.at(id);
            ++id;
        }
    }

    if (id == 0) return 0;

    int c = whichClust.lower_bound(randreal(0, 1))->second;
    vector<int> carr_idx;
    for (int u = 0; u < static_cast<int>(worldData->carriers->at(c).size()); ++u) carr_idx.push_back(u);
    random_shuffle(carr_idx.begin(), carr_idx.end());

    std::shared_ptr<Chromosome> chrom = worldData->carriers->at(c).at(carr_idx.at(0));
    std::shared_ptr<Chromosome> chrom2 = worldData->carriers->at(c).at(carr_idx.at(1));

    if (!(chrom->getContext() == chrom2->getContext())) {
        std::cerr << "Error: attempted coalescence of different contexts in coalesceEvent_SMC. "
                  << "ctx1(pop=" << chrom->getContext().pop << ",inv=" << chrom->getContext().inversion << ") "
                  << "ctx2(pop=" << chrom2->getContext().pop << ",inv=" << chrom2->getContext().inversion << ")\n";
        return 0;
    }

    std::shared_ptr<ARGNode> newNode;
    newNode.reset(new ARGNode(worldData->argNodeVec.size(), chrom, chrom2, worldData->generation));
    worldData->argNodeVec.push_back(newNode);
    chrom->merge(chrom2, newNode);

    worldData->carriers->at(c).erase(remove(worldData->carriers->at(c).begin(),
                                            worldData->carriers->at(c).end(),
                                            chrom2),
                                     worldData->carriers->at(c).end());

    return 1;
}
