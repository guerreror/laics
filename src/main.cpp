/*
 *  InvCoal
 *  Version:  labp_v19 (modified to read YAML params and name gene tree CSV files 
 *          according to the input site position value)
 *
 *  This is the main file for the simulation.
 *  Parameter values are now updated via YAML files using the Python wrapper.
 */

// Includes from STL:
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <chrono>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <numeric>
using namespace std;

#include <memory>
using std::shared_ptr;
using std::unique_ptr;
using std::weak_ptr;

// Includes for our files:
#include "parameters.h"
#include "world.h"
#include "ran_mk.h"
#include "chromosome.h"
#include "sitenode.h"    // SiteNode now includes our new CSV output method.
#include "snptree.h"
#include "tajima.h"
#include "migprob.h"     // <-- NEW: time-scheduled migration matrices

// Global declarations:
std::random_device rd;
auto seed = rd();
std::mt19937_64 gen(seed); // Random generator declared globally.

vector<vector<double>> buildMigMatrix(Parameters &p)
{
    // ======================================================================================
    // A simple stepping-stone migration matrix.
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
                mig_prob.at(i).at(j) = 1 - m; // probability of staying in current pop
            else if (i == j - 1)
            {
                mig_prob.at(i).at(j) = m / 2; // probability of migrating to neighboring pop (to the right)
                if (i == 0)
                    mig_prob.at(i).at(j) = m; // if first pop, adjust probability
            }
            else if (i == j + 1)
            {
                mig_prob.at(i).at(j) = m / 2; // probability of migrating to neighboring pop (to the left)
                if (i == nPops - 1)
                    mig_prob.at(i).at(j) = m; // if last pop, adjust probability
            }
            else
                mig_prob.at(i).at(j) = 0;
        }
    }

    // (Additional mapping code is present in the original, but the returned mig_prob is used.)
    return mig_prob;
}

int main(int argc, const char *argv[])
{
    std::cerr << "Random Seed: " << seed << '\n';

    // If a seed is provided as the first argument, update the seed.
    auto input_seed = std::stoi(argv[1]);
    if (input_seed != -1)
    {
        std::cerr << "Input Seed: " << input_seed << '\n';
        seed = input_seed;
        gen.seed(input_seed);
    }

    // Timer start
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    unsigned int totalEvents = 0;

    stringstream infile, ms_ss, sstat_ss;
    infile << "inLABP.pars";  // (dummy; actual parameters are passed via the Python wrapper)
    ms_ss << "tests/outLABP_" << seed << ".sites";
    sstat_ss << "tests/outLABP_" << seed << ".stats";

    // Collect remaining command-line arguments (starting at argv[2]) as parameter strings.
    vector<string> param_vec(argv + 2, argv + argc);
    // Construct Parameters using these values.
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
        std::cerr << "[mig] loaded " << schedule.size() << " matrices from " << mig_json << "\n";
    } else {
        mig_prob = buildMigMatrix(params);
        std::cerr << "[mig] no JSON; using fallback builder.\n";
    }
    // -----------------------------

    vector<double> outTime(params.paramData->n_SNPs, 0.0);

    ofstream msout, stout;
    if (params.paramData->msOutput)
    {
        msout.open(ms_ss.str().c_str());
        stout.open(sstat_ss.str().c_str());
    }
    if (stout.is_open())
        stout << "ET totL S piTotal pi1 pi2 fst dxy tajD\n";

    double LDsum = 0.0;
    int ticker = (int)(nRuns / 10);
    if (ticker > 0)
        std::cerr << "Progress ticker (1 tick = 10% of runs): ";

    for (int timer = 0; timer < (int)nRuns; ++timer)
    {
        if (ticker > 0 && timer % ticker == 0)
            std::cerr << "+";

        params.setPhi();
        params.setSNPs();
        params.setCarriers();
        unsigned nCarriers = params.paramData->initChr.size();

        World *world = new World(params.getpData());
        while (!world->simulationFinished())
        {
            // -----------------------------
            // NEW: if we have a schedule, switch matrices when generation crosses the next timestamp
            // -----------------------------
            if (!schedule.empty() && next_idx < schedule.size()) {
                double now   = world->nGenerations();
                double tnext = schedule[next_idx].first;
                if (now >= tnext) {
                    mig_prob = schedule[next_idx].second;
                    ++next_idx;
                    std::cerr << "[mig] switched matrix @t=" << tnext
                              << " (rows=" << mig_prob.size() << ")\n";
                }
            }
            // -----------------------------

            totalEvents += world->simulateGeneration(mig_prob);
        }

        vector<shared_ptr<ARGNode>> allNodes = world->getARGVec();
        vector<double> tempLD(params.paramData->n_SNPs, 0.0);
        vector<SNPtree> trees;
        unsigned novar = 0;
        vector<double> varPos;
        double lengthLastSite = 0.0;

        // For each site, construct the gene tree and process SNP statistics.
        for (int k = 0; k < (int)nSites; ++k)
        {
            unsigned pos = k;
            if (!params.paramData->fixedS)
                pos = 0;  // if not fixedS, all sites use the same gene tree

            // Create gene tree for the current site.
            SiteNode geneTree(params.paramData->neut_site[pos], allNodes.back());
            geneTree.calcBranchLengths(0);
            geneTree.calcBranchLengths_informative(0);

            // ---------------------------
            // Write CSV file for gene tree.
            double originalValue = params.paramData->snpPositions[pos] * params.paramData->BasesPerMorgan;
            int bp_val = static_cast<int>(originalValue);
            stringstream csvFileName;
            csvFileName << "genetree_siteposition_" << bp_val << ".csv";
            geneTree.writeCSV(csvFileName.str());
            // ---------------------------

            SNPtree tmp(geneTree, nCarriers, params.paramData->theta);
            unsigned mutantCount = accumulate(tmp.SNPvalues.begin(), tmp.SNPvalues.end(), 0U);
            if (mutantCount > 0)
            {
                trees.push_back(tmp);
                varPos.push_back(k);
            }
            else
            {
                ++novar;
            }
            lengthLastSite = tmp.totalLength / (double)params.paramData->totalPopSize;
            double totalmrca = geneTree.getTime() / (double)params.paramData->totalPopSize;
            outTime.at(k) += totalmrca;
            tempLD[k] = totalmrca;
        }

        double pi = 0.0, piP1 = 0.0, piP2 = 0.0, dxy = 0.0, fst = 0.0;
        double varSites = nSites - novar;
        for (int i = 0; i < (int)varSites; ++i)
        {
            vector<unsigned> snps = trees.at(i).SNPvalues;
            vector<vector<unsigned>> split = split_by_population(snps, params.getSamplePerPop());
            double totalPi = heterozygosity(snps);
            pi += totalPi;
            vector<double> pipop = pi_by_pop(split);
            piP1 += pipop[0];
            piP2 += pipop[1];
            fst += fst_nei(totalPi, pipop);
            dxy += calcdxy(split);
        }
        if (varSites > 0) {
            pi   /= nSites;
            piP1 /= nSites;
            piP2 /= nSites;
            fst  /= varSites;
            dxy  /= nSites;
        }

        vector<double> posout = params.paramData->neut_site;
        if (!params.paramData->fixedS)
            posout = varPos;

        if (msout.is_open()) {
            msout << "// \nsegsites: " << varSites << "\npositions: ";
            for (int k = 0; k < (int)varSites; ++k)
                msout << " " << posout.at(k) << " ";
            msout << "\n";
            for (int i = 0; i < nCarriers; ++i) {
                for (int k = 0; k < (int)varSites; ++k)
                    msout << trees[k].SNPvalues[i];
                msout << "\n";
            }
        }
        double tempProd = tempLD[0] * tempLD[nSites - 1];
        LDsum += tempProd;

        if (stout.is_open()) {
            stout << tempLD[0] << " " << lengthLastSite << " " << varSites << " " 
                  << pi << " " << piP1 << " " << piP2 << " " 
                  << fst << " " << dxy << " " << tajd(nCarriers, varSites, pi) << "\n";
        }

        delete world;
    } // end main simulation loop
    

    std::cerr.precision(6);
    std::cerr << "\nMean LD (E[T1,n]) = " << LDsum / (double)nRuns << "\n";
    if (!params.paramData->fixedS)
        std::cerr << "Mean TMRCA = " << outTime.at(0) / nRuns << "\n";
    else {
        std::cerr << "Mean TMRCA per site (E[T1]...E[Tn])\n";
        for (int k = 0; k < (int)nSites; ++k)
            std::cerr << " " << outTime.at(k) / nRuns;
        std::cerr << "\n";
    }
    if (msout.is_open())
        msout.close();
    if (stout.is_open())
        stout.close();

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cerr << "Elapsed time: " << elapsed_seconds.count() << "s\n";

    return 0;
}
