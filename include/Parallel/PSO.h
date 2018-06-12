/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde and Authors------------
----------------- e-mail: amin.kamal.2016@uni.strath.ac.uk ---------------------
--------------------------- Author: Amin Kamal ---------------------------------

 Original non-OOP C Code written by: Yuhui Shi, yuhui.shi@eds.com, 1998
 http://www.engr.iupui.edu/~eberhart/web/psobook/download/download.html
*/

#ifndef PSO_H
#define PSO_H

#include <vector>
#include "Parallel/EvolutionaryAlgorithms.h"
#include <thread>

namespace Parallel
{
    class PSO : public EvolutionaryAlgorithm
    {
        Population vx;
        Population xx; // xx "distance" is our population, thus the inherited variable population is not used to
                       // maintain the original code's notion
        Population pbestx;
        std::vector<double> pbest;

        int  NUMBER_OF_AGENTS;
        double MAXV, MAXX;
        double weight, weight_up;
        int  MAXITER;
        int DIMENSION;
        int a,b;
        int i, j;
        int iter;
        int gbest;
        double minval;
        bool firsttime;
        long nfeval;

        //For CEC2017
        int func_num;

    public:
        void init();
        void executeEA();
        int getNumberOfGenerations();
        Population getIndividuals(int n, bool remove = true);
        void acceptTake(Population p);
        void acceptFill(Population p);
        void acceptBuffer(Population p);
        void generatePopulationForIMEA();
        void incrementRunCount();
        void output();
        void reset();

        void saveOutputToFile(std::string description, std::string filename);
        Population getPopulation();

        PSO(int d, int numberOfAgents, double maxVelocity, std::vector<std::vector<double>> dbounds,
            int maxIterations, double wt, int cec_func_num);
    };
}



#endif
