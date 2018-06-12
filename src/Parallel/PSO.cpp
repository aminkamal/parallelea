/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde and Authors------------
----------------- e-mail: amin.kamal.2016@uni.strath.ac.uk ---------------------
--------------------------- Author: Amin Kamal ---------------------------------
*/

#include "Parallel/PSO.h"

namespace Parallel
{
    //
    // Implementation of virtual functions from EvolutionaryAlgortihm.h
    // Each EA must implement these
    //

    int PSO::getNumberOfGenerations()
    {
        return iter;
    }

    void PSO::incrementRunCount()
    {
        currentRun++;
    }

    Population PSO::getIndividuals(int n, bool remove)
    {
        pbestx.sort();
        Population tempPopulation;
        tempPopulation.setID(currentlyProcessingSampleID);
        if(remove){
            for(int i=0;i<n;i++)
            {
                Individual ind = xx.front();
                tempPopulation.push_back(ind);
                xx.pop_front();
            }
        }else{
            tempPopulation.insert(tempPopulation.end(), xx.begin(), xx.begin() + n);
        }

        return tempPopulation;
    }

    void PSO::acceptFill(Population p)
    {
        currentlyProcessingSampleID = p.getID();
        xx.insert(xx.end(), p.begin(), p.end());
        init();
    }

    void PSO::acceptTake(Population p)
    {
        currentlyProcessingSampleID = p.getID();
        xx.insert(xx.end(), p.begin(), p.end());
        NUMBER_OF_AGENTS = (int) xx.size();

        init();
    }

    void PSO::acceptBuffer(Population p)
    {
        xx.insert(xx.end(), p.begin(), p.end());
        xx.sort();
        xx.erase(xx.end() - (int)p.size(), xx.end());

        init();
    }

    void PSO::generatePopulationForIMEA()
    {
        xx = EvolutionaryAlgorithm::generatePopulation(DIMENSION, NUMBER_OF_AGENTS, bounds, func_num);
    }

    void PSO::reset()
    {
        xx.clear();
        vx.clear();
        pbest.clear();
        iter = 0;
        nfeval = 0;
        firsttime = true;
    }

    /*
        Start of PSO Code
    */

    PSO::PSO(int d, int numberOfAgents, double maxVelocity, std::vector<std::vector<double>> dbounds,
             int maxIterations, double wt, int cec_func_num)
            : DIMENSION(d), NUMBER_OF_AGENTS(numberOfAgents), MAXV(maxVelocity),
              MAXITER(maxIterations), weight(wt), func_num(cec_func_num),
              nfeval(0), iter(0), minval(0.0), firsttime(true)
    {
        currentRun = 1;
        genmax = maxIterations;
        bounds = dbounds;

        //to avoid each thread initing cec
        if(func_num != 0)
            cec17.cec17_init(d, func_num);
    }

    void PSO::init()
    {
        vx = xx; //fills the deque with the same size as population xx
        pbestx = xx;
        gbest = 0;
        pbest.reserve(NUMBER_OF_AGENTS);
        pbest[gbest] = xx[0].getFitness();

        /* **********************************************
         * This loop initializes the individual agents  for each run
         ********************************************** */
        for (a=0;a<NUMBER_OF_AGENTS;a++)
        {
            for (b=0;b<DIMENSION;b++)
            {
                vx[a][b] = MAXV*(doubleRand(0, 1));
                if (doubleRand(0, 1) > 0.5) vx[a][b]=-vx[a][b];
            }
        }
    }

    void PSO::executeEA()
    {
        if(!seeded){
            seed = (unsigned) time(nullptr) + std::hash<std::thread::id>()(std::this_thread::get_id());
            seeded = true;
        }
        
        /* *******************************************************
            Main Work Loop for each run here
        ******************************************************** */

        iter++;

        //update inertia weight
        weight_up = (weight-0.4) * (MAXITER - iter) /MAXITER +0.4;    //time variant weight, linear from weight to 0.4

        //weight_up=weight;		//constant inertia weight

        for (a=0;a<NUMBER_OF_AGENTS;a++)
        {
            if(func_num != 0)
                cec17.evaluate(xx[a], DIMENSION, 1, func_num);
            else
                evaluateNonCEC(xx[a]);
            nfeval++;
            minval = xx[a].getFitness();

            if (firsttime) pbest[a]=minval;

            if (minval < pbest[a])
            {
                pbest[a]=minval;
                pbestx[a]=xx[a];
                if (pbest[a] < pbest[gbest])
                {
                    gbest=a;
                }
            }

            for (b=0;b<DIMENSION;b++)
            {
                vx[a][b] = weight_up*vx[a][b] + 2*doubleRand(0, 1)*(pbestx[a][b]-xx[a][b]) +
                           2*doubleRand(0, 1)*(pbestx[gbest][b]-xx[a][b]);
                if (vx[a][b]>MAXV)
                    vx[a][b]=MAXV;
                else if (vx[a][b]<-MAXV)
                    vx[a][b]=-MAXV;
            }
            xx[a] = xx[a] + vx[a]; //Define new coordinates
        }

        if(func_num != 0)
            storeCECError(currentRun, iter, func_num, pbest[gbest]);

        firsttime = false;
    }

    void PSO::output()
    {
        printf(
                "PSO\nPSO Parameters:\n"
                        "Number of agents\t%d\n"
                        "Maximum Velocity\t%lf\n"
                        "Maximum Distance\t%lf\n"
                        "Weight factor\t\t%lf\n\n",
                NUMBER_OF_AGENTS, MAXV, MAXX, weight);
        printf("Best cost function value: %f\n\nMembers:\n\n", pbest[gbest]);
        for(int kk = 0; kk < DIMENSION; kk++)
            printf("D%d=%f\n", kk + 1, pbestx[gbest][kk]);
    }

    void PSO::saveOutputToFile(std::string description, std::string filename)
    {
        std::ofstream outfile;
        outfile.open ("output/" + filename + ".txt", std::ios_base::app);
        if(!outfile.good()){
            std::cout << "OUTPUT ERROR: Ensure that a folder named 'output' exists!" << std::endl;
            return;
        }
        outfile << "========================" << std::endl;
        outfile << description << std::endl;
        outfile << "========================" << std::endl;

        outfile << "PSO\nPSO Parameters: \n";
        outfile << "Number of agents\t" << NUMBER_OF_AGENTS << "\n";
        outfile << "Maximum Velocity\t" << MAXV << "\n";
        outfile << "Maximum Distance\t" << MAXX << "\n";
        outfile << "Weight factor\t" << weight << "\n";

        outfile << "Best cost function value: " <<  pbest[gbest] << "\n";
        for(int kk = 0; kk < DIMENSION; kk++)
            outfile << "D" << kk + 1 << "=" << pbestx[gbest][kk] << "\n";

        outfile << "\n\n\n";

        outfile.close();
    }

    Population PSO::getPopulation()
    {
        return pbestx;
    }
}