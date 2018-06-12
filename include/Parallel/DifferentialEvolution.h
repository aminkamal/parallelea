/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde and Authors------------
----------------- e-mail: amin.kamal.2016@uni.strath.ac.uk ---------------------
--------------------------- Author: Amin Kamal ---------------------------------

Original non-OOP DE by Rainer Storn and Ken Price:
http://www1.icsi.berkeley.edu/~storn/code.html
*/

#ifndef DIFFERNTIALEVOLUTION_H
#define DIFFERNTIALEVOLUTION_H

#include <thread>
#include "EvolutionaryAlgorithms.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "memory.h"

namespace Parallel
{
    class DifferentialEvolution : public EvolutionaryAlgorithm
    {
    private:
        Population newpopulation;
        Individual best;
        Individual bestit;
        Individual tmp;

        void crossover();
        void executeEA();

        int   i, j, L, n;      /* counting variables                 */
        int   r1, r2, r3, r4;  /* placeholders for random indexes    */
        int   r5;              /* placeholders for random indexes    */
        int   D;               /* Dimension of parameter vector      */
        double CR;            /* control variables of DE            */
        int   strategy;        /* choice parameter for screen output */
        double F;            /* control variables of DE            */
        int   NP;              /* number of population members       */
        double trial_cost;      /* buffer variable                    */
        long  nfeval;          /* number of function evaluations     */
        int   gen;
        double cvar;            /* computes the cost variance         */
        double cmean;           /* mean cost                          */
        int   imin;            /* index to member with lowest energy */
        double cmin;            /* help variables                     */

        //For CEC2017
        int func_num;

    public:
        void init();
        //For both PoolEA and IMEA
        int getNumberOfGenerations();
        void incrementRunCount();

        //For PoolEA
        Population getIndividuals(int n, bool remove = true);
        void acceptTake(Population p);
        void acceptFill(Population p);
        void acceptBuffer(Population p);

        //For IMEA
        void generatePopulationForIMEA();

        void output();
        void reset();

        void saveOutputToFile(std::string description, std::string filename);
        Population getPopulation();

        DifferentialEvolution(int d, int gmax, double cr, int s, double f, int np,
                              std::vector<std::vector<double>> dbounds, int cec_func_num);

    };
}
#endif
