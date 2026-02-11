/*
 *  smc.cpp
 *  SMC entry point (minimal stub).
 *
 *  Runs the standard setup and ARG simulation, then emits the first
 *  site tree and exits. 
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <chrono>
#include <random>

using namespace std;

#include "parameters.h"
#include "world.h"
#include "ran_mk.h"
#include "sitenode.h"
#include "snptree.h"
#include "migprob.h"

std::random_device rd;
auto seed = rd();
std::mt19937_64 gen(seed);

static vector<vector<double>> buildMigMatrix(Parameters &p)
{
    unsigned int nPops = static_cast<int>(p.paramData->popSizeVec.size());
    vector<vector<double>> mig_prob;
    mig_prob.resize(nPops);
    double m = p.paramData->migRate.at(0) / (2 * p.paramData->totalPopSize);

    for (unsigned int i = 0; i < nPops; i++)
        mig_prob.at(i).resize(nPops);
    for (unsigned int i = 0; i < nPops; i++)
    {
        for (unsigned int j = 0; j < nPops; j++)
        {
            if (nPops == 1)
                mig_prob.at(i).at(j) = 1;
            else if (i == j)
                mig_prob.at(i).at(j) = 1 - m;
            else if (i == j - 1)
            {
                mig_prob.at(i).at(j) = m / 2;
                if (i == 0)
                    mig_prob.at(i).at(j) = m;
            }
            else if (i == j + 1)
            {
                mig_prob.at(i).at(j) = m / 2;
                if (i == nPops - 1)
                    mig_prob.at(i).at(j) = m;
            }
            else
                mig_prob.at(i).at(j) = 0;
        }
    }
    return mig_prob;
}

int main(int argc, const char *argv[])
{
    std::cerr << "Random Seed: " << seed << '\n';

    auto input_seed = std::stoi(argv[1]);
    if (input_seed != -1)
    {
        std::cerr << "Input Seed: " << input_seed << '\n';
        seed = input_seed;
        gen.seed(input_seed);
    }

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    stringstream infile;
    infile << "inLABP.pars";

    vector<string> param_vec(argv + 2, argv + argc);
    Parameters params(infile.str().c_str(), param_vec);

    unsigned int nRuns = params.paramData->nRuns;
    unsigned int nSites = params.paramData->n_SNPs;

    const std::string mig_json = "src/migration_matrices.json";
    auto schedule = readMigrationSchedule(mig_json);

    vector<vector<double>> mig_prob;
    size_t next_idx = 0;

    if (!schedule.empty()) {
        double g0 = 0.0;
        while (next_idx < schedule.size() && schedule[next_idx].first <= g0) {
            mig_prob = schedule[next_idx].second;
            ++next_idx;
        }
        if (mig_prob.empty()) {
            mig_prob = schedule.front().second;
        }
    } else {
        mig_prob = buildMigMatrix(params);
    }

    for (int timer = 0; timer < (int)nRuns; ++timer)
    {
        params.setPhi();
        params.setSNPs();
        params.setCarriers();

        World *world = new World(params.getpData());
        while (!world->simulationFinished())
        {
            if (!schedule.empty() && next_idx < schedule.size()) {
                double now   = world->nGenerations();
                double tnext = schedule[next_idx].first;
                if (now >= tnext) {
                    mig_prob = schedule[next_idx].second;
                    ++next_idx;
                }
            }
            world->simulateGeneration(mig_prob);
        }

        vector<shared_ptr<ARGNode>> allNodes = world->getARGVec();
        if (nSites == 0 || allNodes.empty()) {
            std::cerr << "No sites or ARG nodes available.\n";
            delete world;
            break;
        }

        unsigned pos = 0;
        if (!params.paramData->fixedS)
            pos = 0;

        SiteNode geneTree(params.paramData->neut_site[pos], allNodes.back());
        geneTree.calcBranchLengths(0);
        geneTree.calcBranchLengths_informative(0);
        geneTree.writeCSV("genetree_first_site.csv");
        geneTree.writeDOT("genetree_first_site.dot");

        double std_len = 0.0;
        double inv_len = 0.0;
        geneTree.getTotalLengthByInversion(std_len, inv_len);
        std::cerr << "Ls(x) = " << std_len << "\n";
        std::cerr << "Li(x) = " << inv_len << "\n";
        std::cerr << "L(x) = " << (std_len + inv_len) << "\n";
        std::cerr << "L(x) (getTotalLength func) = "
                  << geneTree.getTotalLength(0.0) << "\n";

        double r = 1.0 / params.paramData->BasesPerMorgan;
        vector<EdgeWeight> standard_edges;
        vector<EdgeWeight> inverted_edges;
        geneTree.getEdgeWeightsByInversion(params.paramData->popSizeVec,
                                           params.paramData->initialFreqs,
                                           r,
                                           standard_edges,
                                           inverted_edges);

        {
            std::ofstream std_out("edge_weights_standard.csv");
            if (std_out.is_open()) {
                std_out << "parent,child,weight\n";
                for (const auto &e : standard_edges) {
                    std_out << e.parent << "," << e.child << "," << e.weight << "\n";
                }
            }
        }
        {
            std::ofstream inv_out("edge_weights_inverted.csv");
            if (inv_out.is_open()) {
                inv_out << "parent,child,weight\n";
                for (const auto &e : inverted_edges) {
                    inv_out << e.parent << "," << e.child << "," << e.weight << "\n";
                }
            }
        }

        delete world;
        break;
    }

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cerr << "Elapsed time: " << elapsed_seconds.count() << "s\n";

    return 0;
}
