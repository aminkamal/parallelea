/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde and Authors------------
----------------- e-mail: amin.kamal.2016@uni.strath.ac.uk ---------------------
--------------------------- Author: Amin Kamal ---------------------------------
*/

#ifndef EVOLUTIONARYALGORITHMS_H
#define EVOLUTIONARYALGORITHMS_H

#include <string>
#include <map>
#include <deque>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "Parallel/Population.h"
#include "Parallel/cec17_test_func.h"

namespace Parallel
{
    enum class ExecutionMethod {Sync, Async};
    class EvolutionaryAlgorithm
    {
    private:
        static CEC2017Problems scec17;
    protected:
        Population population;
        std::vector<std::vector<double>> bounds;
        CEC2017Problems cec17;
        unsigned int seed;
        bool seeded = false;

        int currentlyProcessingSampleID;    //For PoolEA
        int currentRun;                     //For CEC2017 problems
        int genmax;
        void storeCECError(int run, int gen, int func_num, double cmin);
        static const int multipliersCount = 14;
        double multipliers[multipliersCount] = {0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
        double doubleRand(double min, double max) {
            return (min + (double) rand_r(&seed)/RAND_MAX * (max - min));
        }
    public:
        virtual void init() = 0;
        virtual void executeEA() = 0;
        virtual int getNumberOfGenerations() = 0;
        virtual Population getIndividuals(int n, bool remove = true) = 0;
        virtual void acceptTake(Population p) = 0;
        virtual void acceptFill(Population p) = 0;
        virtual void acceptBuffer(Population p) = 0;
        virtual void generatePopulationForIMEA() = 0;
        virtual void incrementRunCount() = 0;
        virtual void output() = 0;
        virtual void reset() = 0;
        virtual void saveOutputToFile(std::string description, std::string filename) = 0;
        virtual Population getPopulation() = 0;
        //CEC2017
        double CEC2017ErrorValues[100][14] = {};
        void printCECError();
        void saveCECErrorToFile(std::string description, std::string filename);
        static unsigned int s_seed;
        static double s_doubleRand(double min, double max) {
            return (min + (double) rand_r(&EvolutionaryAlgorithm::s_seed)/RAND_MAX * (max - min));
        }

        static Population generatePopulation(int D, int populationSize, std::vector<std::vector<double>> bounds,
                                             int func_num);

        virtual ~EvolutionaryAlgorithm() {};
    };

    extern void evaluateNonCEC(Individual &ind);
}

#endif
