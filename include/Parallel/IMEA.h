/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde and Authors------------
----------------- e-mail: amin.kamal.2016@uni.strath.ac.uk ---------------------
--------------------------- Author: Amin Kamal ---------------------------------
*/

#ifndef IMEA_H
#define IMEA_H

#include "Graph.h"
#include <vector>
#include <memory>
#include "EvolutionaryAlgorithms.h"
#include <thread>
#include <mutex>
#include <atomic>

namespace Parallel
{
    class IMEA;

    enum class Topology{
        Ring, FullyConnected, Broadcast
    };

    class Island
    {
    private:
        int islandID_;
        IMEA &imea_;
        Population buffer;

    public:
        EvolutionaryAlgorithm & ea_;
        Island(int islandID, IMEA &imea, EvolutionaryAlgorithm &ea) : islandID_(islandID), imea_(imea), ea_(ea) {}
        int getID();
        bool isFinished();
        void initEA();
        void appendToBuffer(Population p);
        void processBuffer();

        void executeSync();
        void executeAsync();
    };

    class IMEA
    {
    private:
        UndirectedGraph g;
        int demeSize_;
        int workerGenerations_;
        int numberOfRuns_;
        int func_num_;
        int migrationRate_;
        std::vector<Island> islV;
        const int maxNumberOfThreads = 32;

    public:
        IMEA(int demeSize, int workerGenerations, int numberOfRuns, int migrationRate, int func_num = 0)
                : demeSize_(demeSize),
                  workerGenerations_(workerGenerations),
                  numberOfRuns_(numberOfRuns),
                  func_num_(func_num),
                  migrationRate_(migrationRate) {}
        void addIsland(Island isl);
        void start(ExecutionMethod em, std::string filename, bool printError = false);

        void setTopology(Topology t, int numberOfIslands);
        void spreadPopulation(int islandID, Population p);

        int getMaxNumberOfGenerations();
        int getMigrationRate();
        int getDemeSize();
        //Threading
        std::mutex islandMutex;
        std::atomic<int> syncFinishCounter;
    };

}

#endif
