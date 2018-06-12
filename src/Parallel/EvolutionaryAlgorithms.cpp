/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde and Authors------------
----------------- e-mail: amin.kamal.2016@uni.strath.ac.uk ---------------------
--------------------------- Author: Amin Kamal ---------------------------------
*/

#include "Parallel/EvolutionaryAlgorithms.h"

namespace Parallel
{
    CEC2017Problems EvolutionaryAlgorithm::scec17;
    unsigned int EvolutionaryAlgorithm::s_seed;

    Population EvolutionaryAlgorithm::generatePopulation(int D, int populationSize,
                                                         std::vector<std::vector<double>> bounds, int func_num)
    {
        s_seed = (unsigned) time(nullptr);
        //to avoid each thread initing cec
        if(func_num != 0)
            scec17.cec17_init(D, func_num);

        Population temp;
        for(int i = 0;i<populationSize;i++){
            std::vector<double> members;
            double fitness = 0.0;
            for(int j = 0; j < D; j++)
            {
                if(bounds.size() == 1)
                    members.push_back(s_doubleRand(bounds[0][0], bounds[0][1])); //0 = Lower bound, 1 = Upper bound
                else
                    members.push_back(s_doubleRand(bounds[j][0], bounds[j][1]));
            }
            Individual ind(i, members, fitness);
            if(func_num != 0)
                scec17.evaluate(ind, D, 1, func_num);
            else
                evaluateNonCEC(ind);
            temp.push_back(ind);
        }
        return temp;
    }

    void EvolutionaryAlgorithm::printCECError()
    {
        for(int i = 0; i < currentRun - 1; i++) //currentRun is at max after the algorithm has finished
            std::cout << "RUN" << i + 1
                      << ", ERROR" << i + 1 << "="
                      << CEC2017ErrorValues[i][multipliersCount - 1]
                      << std::endl;
    }

    void EvolutionaryAlgorithm::saveCECErrorToFile(std::string description, std::string filename)
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
        for(int i = 0; i < currentRun - 1; i++)
            for(int j = 0; j < multipliersCount; j++)
                outfile << "RUN" << i + 1 << ", "
                        << "ERROR" << j + 1 << "\t(Recorded after " << multipliers[j] * genmax << "/" << genmax << ") generations \t="
                        << CEC2017ErrorValues[i][j] << std::endl;
        outfile.close();
    }

    void EvolutionaryAlgorithm::storeCECError(int run, int gen, int func_num, double cmin)
    {
        for(int i=0;i<multipliersCount;i++)
            if(gen >= genmax * multipliers[i] && CEC2017ErrorValues[run-1][i] == 0)
                CEC2017ErrorValues[run-1][i] = cmin - (func_num*100);
    }
}