/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde and Authors------------
----------------- e-mail: amin.kamal.2016@uni.strath.ac.uk ---------------------
--------------------------- Author: Amin Kamal ---------------------------------
*/

#ifndef POOLEA_H
#define POOLEA_H

#include <string>
#include <queue>
#include <deque>
#include <vector>
#include <map>
#include <memory>
#include "EvolutionaryAlgorithms.h"
#include "DifferentialEvolution.h"
#include <thread>
#include <mutex>
#include <atomic>

namespace Parallel
{
    class EvoSpaceCPP; //yes, we have cyclic dependency

    enum class RESPAWNTYPE{
        REINSERT, RANDOM
    };

    enum class TAKETYPE{
        RANDOM, BEST
    };

    class EvoSpaceWorker
    {
    public:
        EvolutionaryAlgorithm & ea_;
        EvoSpaceCPP &espp_;
        EvoSpaceWorker(EvolutionaryAlgorithm &ea, EvoSpaceCPP &espp)
                : ea_(ea), espp_(espp) {}
        bool isFinished();
        void take(int n);
        void storePopulationIntoPool(int n);

        void executeSync();
        void executeAsync();
    };

    class EvoSpaceCPP
    {
    private:
        std::deque<Population> sample_queue;
        std::vector<EvoSpaceWorker> eswV;
        const int maxNumberOfThreads = 32;

        int populationSize_;
        int D_;
        std::vector<std::vector<double>> bounds;
        int sampleSize_;
        int workerGenerations_;
        int numberOfRuns_;
        int migrationRate_;

        int func_num_;
        TAKETYPE TAKET;

        unsigned int MIN_SIZE;
        int RE_INSERT_SAMPLES;
        RESPAWNTYPE RESPAWN;

        int sample_counter;
        int returned_counter;

    public:
        Population pool_population;
        EvoSpaceCPP(int populationSize,
                    int D,
                    std::vector<std::vector<double>> dbounds,
                    int sampleSize,
                    int workerGenerations,
                    int numberOfRuns,
                    int migrationRate,
                    int func_num = 0,
                    TAKETYPE takeType = TAKETYPE ::BEST,
                    unsigned int minSize = 30,
                    int reinsertSamples = 8,
                    RESPAWNTYPE respawnType = RESPAWNTYPE::REINSERT
                    )
                : populationSize_(populationSize),
                  D_(D),
                  sampleSize_(sampleSize),
                  workerGenerations_(workerGenerations),
                  numberOfRuns_(numberOfRuns),
                  migrationRate_(migrationRate),
                  func_num_(func_num),
                  TAKET(takeType),
                  MIN_SIZE(minSize), //Different naming schemes below to keep consistent with the original EvoSpace.py
                  RE_INSERT_SAMPLES(reinsertSamples),
                  RESPAWN(respawnType),
                  sample_counter(0),
                  returned_counter(0)
        {
            bounds = dbounds;
            pool_population = EvolutionaryAlgorithm::generatePopulation(D, populationSize, bounds, func_num);
        }

        int getMaxNumberOfGenerations();
        int getMigrationRate();
        int getSampleSize();
        Population get_sample(int size);
        void put_sample(Population population);
        void respawn_sample(Population population);
        void respawn(int n = 1);
        void addWorker(EvoSpaceWorker esw);
        void start(ExecutionMethod em, std::string filename, bool printError = false);
        void distributeFairly(int n, std::vector<EvoSpaceWorker>& eswV);
        //Threading
        std::mutex poolMutex;
        std::atomic<int> syncFinishCounter;
    };
}

#endif
