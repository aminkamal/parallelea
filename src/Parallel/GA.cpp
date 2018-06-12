/*

                  This file is licensed under GNU LGPL v3
 You should have received a copy of the GNU Lesser General Public License
 along with GA.cpp and GA.h.  If not, see <http://www.gnu.org/licenses/>

------------Copyright (C) 2017 University of Strathclyde and Authors------------
----------------- e-mail: amin.kamal.2016@uni.strath.ac.uk ---------------------
--------------------------- Author: Amin Kamal ---------------------------------
*/

#include "Parallel/GA.h"

namespace Parallel
{
    //
    // Implementation of virtual functions from EvolutionaryAlgortihm.h
    // Each EA must implement these
    //

    int GA::getNumberOfGenerations()
    {
        return iter;
    }

    void GA::incrementRunCount()
    {
        currentRun++;
    }

    Population GA::getIndividuals(int n, bool remove)
    {
        population.sort();
        Population tempPopulation;
        tempPopulation.setID(currentlyProcessingSampleID);
        if(remove){
            for(int i=0;i<n;i++)
            {
                Individual ind = population.front();
                ind.setFitness(ind.getFitness() *-1); //GA-specific
                tempPopulation.push_back(ind);
                population.pop_front();
            }
        }else{
            for(int i=0;i<n;i++)
            {
                Individual ind = population[i];
                ind.setFitness(ind.getFitness() *-1);
                tempPopulation.push_back(ind);
            }
        }
        return tempPopulation;
    }

    void GA::acceptFill(Population p)
    {
        //reverse the sign as this alg maximizes by default
        for(auto &pp : p)
            pp.setFitness(pp.getFitness() *-1);

        currentlyProcessingSampleID = p.getID();
        population.insert(population.end(), p.begin(), p.end());
    }

    void GA::acceptTake(Population p)
    {
        for(auto &pp : p)
            pp.setFitness(pp.getFitness() *-1);

        currentlyProcessingSampleID = p.getID();
        population.insert(population.end(), p.begin(), p.end());
        NP =(int) population.size();
    }

    void GA::acceptBuffer(Population p)
    {
        for(auto &pp : p)
            pp.setFitness(pp.getFitness() *-1);

        population.insert(population.end(), p.begin(), p.end());
        population.sort();
        population.erase(population.end() -(int)p.size(), population.end());
    }

    void GA::generatePopulationForIMEA()
    {
        population = EvolutionaryAlgorithm::generatePopulation(D_, NP, bounds, func_num);
        for(auto &p : population)
            p.setFitness(p.getFitness() *-1);
    }

    void GA::reset()
    {
        population.clear();
        iter = 0;
        nfeval = 0;
    }

    /*
        Start of GA Code
    */

    GA::GA(int d, int populationSize, double CR, double mutationRatio, int maxIterations,
           std::vector<std::vector<double>> dbounds, int cec_func_num)
            : D_(d), NP(populationSize), crossOverRatio_(CR), mutationRatio_(mutationRatio),
              func_num(cec_func_num), nfeval(0), iter(0), bestind(0, {}, 0.0)
    {
        currentRun = 1;
        genmax = maxIterations;
        bounds = dbounds;

        //to avoid each thread initing cec
        if(func_num != 0)
            cec17.cec17_init(d, func_num);
    }

    void GA::init()
    {

    }

    void GA::executeEA()
    {
        if(!seeded){
            seed = (unsigned) time(nullptr) + std::hash<std::thread::id>()(std::this_thread::get_id());
            seeded = true;
        }

        iter++;
        bestind = population[0];

        keep_the_best();

        population.sort(); //replaces selector();
        crossover();
        mutate();
        evaluatePopulation();
        elitist();
        if(func_num != 0)
            storeCECError(currentRun, iter, func_num, bestind.getFitness()*-1);
    }

    void GA::output()
    {
        printf(
                "GA\nGA Parameters:\n"
                        "Population size\t%d\n"
                        "Crossover ratio\t%lf\n"
                        "Mutation ratio\t%lf\n\n",
                NP, crossOverRatio_, mutationRatio_);

        std::cout << "\n";
        std::cout << "  Best member after " << genmax << " generations:\n";
        std::cout << "\n";

        for(i = 0; i < D_; i++)
        {
            std::cout << "  var(" << i << ") = " << bestind[i] << "\n";
        }

        std::cout << "\n";
        std::cout << "  Best fitness = " << bestind.getFitness() * -1 << "\n";
    }

    void GA::saveOutputToFile(std::string description, std::string filename)
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

        outfile << "GA\nGA Parameters: \n";
        outfile << "Population size\t" << NP << "\n";
        outfile << "Crossover ratio\t" << crossOverRatio_ << "\n";
        outfile << "Mutation ratio\t" << mutationRatio_ << "\n";

        outfile << "\n";
        outfile << "  Best member after " << genmax << " generations:\n";
        outfile << "\n";

        for(i = 0; i < D_; i++)
        {
            outfile << "  var(" << i << ") = " << bestind[i] << "\n";
        }

        outfile << "\n";
        outfile << "  Best fitness = " << bestind.getFitness() * -1 << "\n";

        outfile << "\n\n\n";

        outfile.close();
    }

    Population GA::getPopulation()
    {
        Population popcopy = population;
        for(auto &pp : popcopy)
            pp.setFitness(pp.getFitness() *-1);
        return popcopy;
    }

    void GA::evaluatePopulation()
    {
        for(int member = 0; member < NP; member++)
        {
            if(func_num != 0)
                cec17.evaluate(population[member], D_, 1, func_num);
            else
                evaluateNonCEC(population[member]);
            nfeval++;
            population[member].setFitness(population[member].getFitness()*-1); //Remember, this alg maximizes by default
        }
    }

//****************************************************************************80

    void GA::crossover()

//****************************************************************************80
//
//  Purpose:
//
//    CROSSOVER selects two parents for the single point crossover.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, int FIRST, is a count of the number of members chosen.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
    {
        const double a = 0.0;
        const double b = 1.0;
        int mem;
        int one;
        int first = 0;
        double x;

        for(mem = 0; mem < NP; ++mem)
        {
            x = doubleRand(a, b);

            if(x < crossOverRatio_)
            {
                ++first;

                if(first % 2 == 0)
                {
                    Xover(one, mem);
                }
                else
                {
                    one = mem;
                }

            }
        }
    }
//****************************************************************************80

    void GA::elitist()

//****************************************************************************80
//
//  Purpose:
//
//    ELITIST stores the best member of the previous generation.
//
//  Discussion:
//
//    The best member of the previous generation is stored as
//    the last in the array. If the best member of the current
//    generation is worse then the best member of the previous
//    generation, the latter one would replace the worst member
//    of the current population.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, double BEST, the best fitness value.
//
//    Local, double WORST, the worst fitness value.
//
    {
        int i;
        double best;
        int best_mem;
        double worst;
        int worst_mem;

        best = population[0].getFitness();
        worst = population[0].getFitness();

        for(i = 0; i < NP - 1; ++i)
        {
            if(population[i+1].getFitness() < population[i].getFitness())
            {
                if(best <= population[i].getFitness())
                {
                    best = population[i].getFitness();
                    best_mem = i;
                }
                if(population[i+1].getFitness() <= worst)
                {
                    worst = population[i+1].getFitness();
                    worst_mem = i + 1;
                }
            }
            else
            {
                if(population[i].getFitness() <= worst)
                {
                    worst = population[i].getFitness();
                    worst_mem = i;
                }
                if(best <= population[i+1].getFitness())
                {
                    best = population[i+1].getFitness();
                    best_mem = i + 1;
                }
            }

        }
//
//  If the best individual from the new population is better than
//  the best individual from the previous population, then
//  copy the best from the new population; else replace the
//  worst individual from the current population with the
//  best one from the previous generation
//
        if(bestind.getFitness() <= best)
        {
            bestind = population[best_mem];
            bestind.setFitness(population[best_mem].getFitness());
        }
        else
        {
            population[worst_mem] = bestind;
            population[worst_mem].setFitness(bestind.getFitness());
        }
    }

//****************************************************************************80

    void GA::keep_the_best()

//****************************************************************************80
//
//  Purpose:
//
//    KEEP_THE_BEST keeps track of the best member of the population.
//
//  Discussion:
//
//    Note that the last entry in the array Population holds a
//    copy of the best individual.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 December 2007
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, int CUR_BEST, the index of the best individual.
//
    {
        int cur_best = 0;

        for(int mem = 0; mem < NP; mem++)
        {
            if(bestind.getFitness() < population[mem].getFitness())
            {
                cur_best = mem;
                bestind.setFitness(population[mem].getFitness());
            }
        }
//
//  Once the best member in the population is found, copy the genes.
//
        bestind = population[cur_best];
    }
//****************************************************************************80

    void GA::mutate()

//****************************************************************************80
//
//  Purpose:
//
//    MUTATE performs a random uniform mutation.
//
//  Discussion:
//
//    A variable selected for mutation is replaced by a random value
//    between the lower and upper bounds of this variable.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number generator.
//
    {
        const double a = 0.0;
        const double b = 1.0;
        int i;
        int j;
        double x;

        for(i = 0; i < NP; i++)
        {
            for(j = 0; j < D_; j++)
            {
                x = doubleRand(a, b);
                if(x < mutationRatio_)
                {
                    if(bounds.size() == 1)
                        population[i][j] = doubleRand(bounds[0][0], bounds[0][1]);
                    else
                        population[i][j] = doubleRand(bounds[j][0], bounds[j][1]);
                }
            }
        }
    }

    void GA::Xover(int one, int two)

//****************************************************************************80
//
//  Purpose:
//
//    XOVER performs crossover of the two selected parents.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 April 2014
//
//  Author:
//
//    Original version by Dennis Cormier and Sita Raghavan.
//    This C++ version by John Burkardt.
//
//  Local parameters:
//
//    Local, int point, the crossover point.
//
//  Parameters:
//
//    Input, int ONE, TWO, the indices of the two parents.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
    {
        int i;
        int point;
        double t;
//
//  Select the crossover point.
//
        point = (int)doubleRand(0, D_ - 1);
//
//  Swap genes in positions 0 through POINT-1.
//
        for(i = 0; i < point; i++)
        {
            t                  = population[one][i];
            population[one][i] = population[two][i];
            population[two][i] = t;
        }
    }
}