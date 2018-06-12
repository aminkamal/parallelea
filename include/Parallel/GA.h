/*

                  This file is licensed under GNU LGPL v3
 You should have received a copy of the GNU Lesser General Public License
 along with GA.cpp and GA.h.  If not, see <http://www.gnu.org/licenses/>

------------Copyright (C) 2017 University of Strathclyde and Authors------------
----------------- e-mail: amin.kamal.2016@uni.strath.ac.uk ---------------------
--------------------------- Author: Amin Kamal ---------------------------------

 SIMPLE_GA
 Original Code written by: John Burkardt
 https://people.sc.fsu.edu/~jburkardt/cpp_src/simple_ga/simple_ga.html
*/

#ifndef GA_H
#define GA_H

#include "Parallel/EvolutionaryAlgorithms.h"
#include <thread>

#include <iostream>
#include <iomanip>

namespace Parallel
{
    class GA : public EvolutionaryAlgorithm
    {
        Individual bestind;

        int NP;
        int D_;
        int i;
        double crossOverRatio_;
        double mutationRatio_;
        long nfeval;
        int iter;

        //For CEC2017
        int func_num;

    public:
        void crossover();
        void elitist();
        void evaluatePopulation();
        void keep_the_best();
        void mutate();
        void Xover(int one, int two);

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

        GA(int d, int populationSize, double CR, double mutationRatio, int maxIterations,
               std::vector<std::vector<double>> dbounds, int cec_func_num);
    };
}



#endif
