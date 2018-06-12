/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde and Authors------------
----------------- e-mail: amin.kamal.2016@uni.strath.ac.uk ---------------------
--------------------------- Author: Amin Kamal ---------------------------------
*/

#include "Parallel/DifferentialEvolution.h"

namespace Parallel
{
    //
    // Implementation of virtual functions from EvolutionaryAlgortihm.h
    // Each EA must implement these
    //

    int DifferentialEvolution::getNumberOfGenerations()
    {
        return gen;
    }

    void DifferentialEvolution::incrementRunCount()
    {
        currentRun++;
    }

    Population DifferentialEvolution::getIndividuals(int n, bool remove)
    {
        population.sort();
        Population tempPopulation;
        tempPopulation.setID(currentlyProcessingSampleID);
        if(remove){
            for(int i=0;i<n;i++)
            {
                Individual ind = population.front();
                tempPopulation.push_back(ind);
                population.pop_front();
            }
        }else{
            tempPopulation.insert(tempPopulation.end(), population.begin(), population.begin() + n);
        }

        return tempPopulation;
    }

    void DifferentialEvolution::acceptFill(Population p)
    {
        currentlyProcessingSampleID = p.getID();
        population.insert(population.end(), p.begin(), p.end());
        init();
    }

    void DifferentialEvolution::acceptTake(Population p)
    {
        currentlyProcessingSampleID = p.getID();
        population.insert(population.end(), p.begin(), p.end());
        NP = (int) population.size();

        init();
    }

    void DifferentialEvolution::acceptBuffer(Population p)
    {
        population.insert(population.end(), p.begin(), p.end());
        population.sort();
        population.erase(population.end() - (int)p.size(), population.end());

        init();
    }

    void DifferentialEvolution::generatePopulationForIMEA()
    {
        population = EvolutionaryAlgorithm::generatePopulation(D, NP, bounds, func_num);
    }

    void DifferentialEvolution::reset()
    {
        nfeval = 0;
        gen = 0;
        population.clear();
        newpopulation.clear();
    }

    /*
        Start of DE Code
    */

    //
    // Used by both:
    // D            = number of parameters (problem dimension)
    // genmax       = maximum number of generations
    // CR           = crossing over factor
    // strategy     = choice of strategy
    // F            = weight factor
    //
    // Used by IMEA only:
    // NP               = population size
    // inibound_l (lb)  = upper parameter bound for init
    // inibound_h (hb)  = lower parameter bound for init

    DifferentialEvolution::DifferentialEvolution(int d, int gmax, double cr, int s, double f, int np,
                                                 std::vector<std::vector<double>> dbounds, int cec_func_num)
            : D(d), CR(cr), strategy(s),
              F(f), NP(np), nfeval(0), gen(0),
              func_num(cec_func_num), best(0, {}, 0.0), bestit(0, {}, 0.0), tmp(0, {}, 0.0)
    {
        currentRun = 1;
        genmax = gmax;
        bounds = dbounds;

        //to avoid each thread initing cec
        if(func_num != 0)
            cec17.cec17_init(d, func_num);
    }

    void DifferentialEvolution::init()
    {
        cmin = population[0].getFitness();
        imin = 0;
        for (i=1; i<NP; i++)
        {
            if (population[i].getFitness()<cmin)
            {
                cmin = population[i].getFitness();
                imin = i;
            }
        }

        best = population[imin];            /* save best member ever          */
        bestit = population[imin];          /* save best member of generation */

        newpopulation = population;
    }

    void DifferentialEvolution::executeEA()
    {
        if(!seeded){
            seed = (unsigned) time(nullptr) + std::hash<std::thread::id>()(std::this_thread::get_id());
            seeded = true;
        }
        gen++;
        imin = 0;

        for (i=0; i<NP; i++)         /* Start of loop through ensemble  */
        {
            do                        /* Pick a random population member */
            {                         /* Endless loop for NP < 2 !!!     */
                r1 = (int)(doubleRand(0,1)*NP);
            }while(r1==i);

            do                        /* Pick a random population member */
            {                         /* Endless loop for NP < 3 !!!     */
                r2 = (int)(doubleRand(0,1)*NP);
            }while((r2==i) || (r2==r1));

            do                        /* Pick a random population member */
            {                         /* Endless loop for NP < 4 !!!     */
                r3 = (int)(doubleRand(0,1)*NP);
            }while((r3==i) || (r3==r1) || (r3==r2));

            do                        /* Pick a random population member */
            {                         /* Endless loop for NP < 5 !!!     */
                r4 = (int)(doubleRand(0,1)*NP);
            }while((r4==i) || (r4==r1) || (r4==r2) || (r4==r3));

            do                        /* Pick a random population member */
            {                         /* Endless loop for NP < 6 !!!     */
                r5 = (int)(doubleRand(0,1)*NP);
            }while((r5==i) || (r5==r1) || (r5==r2) || (r5==r3) || (r5==r4));


            crossover();

            /*=======Trial mutation now in tmp[]. Test how good this choice really was.==================*/

            /* Evaluate new vector in tmp[] */
            if(func_num != 0)
                cec17.evaluate(tmp, D, 1, func_num);
            else
                evaluateNonCEC(tmp);
            nfeval++;
            trial_cost = tmp.getFitness();

            if (trial_cost <= population[i].getFitness())   /* improved objective function value ? */
            {
                population[i].setFitness(trial_cost);
                newpopulation[i] = tmp;
                if (trial_cost<cmin)          /* Was this a new minimum? */
                {                               /* if so...*/
                    cmin=trial_cost;           /* reset cmin to new low...*/
                    imin=i;
                    best = tmp;
                }
            }
            else
            {
                newpopulation[i] = population[i]; /* replace target with old value */
            }
        }   /* End mutation loop through pop. */

        bestit = best;  /* Save best population member of current iteration */

        /* swap population arrays. New generation becomes old one */
        newpopulation.swap(population);

        if(func_num != 0)
            storeCECError(currentRun, gen, func_num, cmin);
    }

    void DifferentialEvolution::output()
    {
        printf(
                "Differential Evolution\nDE Parameters:\n"
                "Population size\t%d\n"
                "Crossover ratio\t%lf\n"
                "Weight factor\t%lf\n"
                "Strategy\t%d\n",
                NP, CR, F, strategy);
        /*----Compute the energy variance (just for monitoring purposes)-----------*/
        cmean = 0.;          /* compute the mean value first */
        for (j=0; j<NP; j++)
        {
            cmean += population[j].getFitness();
        }
        cmean = cmean/NP;

        cvar = 0.;           /* now the variance              */
        for (j=0; j<NP; j++)
        {
            cvar += (population[j].getFitness() - cmean)*(population[j].getFitness() - cmean);
        }
        cvar = cvar/(NP-1);

        printf("\n\nDifferential Evolution\n Best cost funct. value=%-15.10g\n",cmin);
        for (j=0;j<D;j++)
        {
            printf("\n best[%d]=%-15.10g",j,best[j]);
        }
        printf("\n\n Generation=%d  NFEs=%ld   Strategy: %1d    ",gen,nfeval,strategy);
        printf("\n NP=%d    F=%-4.2g    CR=%-4.2g   cost-variance=%-10.5g\n\n", NP,F,CR,cvar);
    }

    void DifferentialEvolution::saveOutputToFile(std::string description, std::string filename)
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

        outfile << "Differential Evolution\nDE Parameters:\n";
        outfile << "Population size\t" << NP << "\n";
        outfile << "Crossover ratio\t" << CR << "\n";
        outfile << "Weight factor\t" << F << "\n";
        outfile << "Strategy\t" << strategy << "\n";

        /*----Compute the energy variance (just for monitoring purposes)-----------*/
        cmean = 0.;          /* compute the mean value first */
        for (j=0; j<NP; j++)
        {
            cmean += population[j].getFitness();
        }
        cmean = cmean/NP;

        cvar = 0.;           /* now the variance              */
        for (j=0; j<NP; j++)
        {
            cvar += (population[j].getFitness() - cmean)*(population[j].getFitness() - cmean);
        }
        cvar = cvar/(NP-1);

        outfile << "\nBest cost funct. value=" << cmin;
        for (j=0;j<D;j++)
        {
            outfile << "\nbest[" << j << "]=" << best[j];
        }
        outfile << "\n Generation=" << gen << "    NFEs=" << nfeval << "    Strategy: " << strategy;
        outfile << "\n NP=" << NP << "    F=" << F << "    CR= " << CR << "    cost-variance=" << cvar;

        outfile << "\n\n\n";

        outfile.close();
    }

    Population DifferentialEvolution::getPopulation()
    {
        return population;
    }

    void DifferentialEvolution::crossover()
    {
        /*=======EXPONENTIAL CROSSOVER============================================================*/

        /*-------DE/best/1/exp--------------------------------------------------------------------*/
        /*-------Our oldest strategy but still not bad. However, we have found several------------*/
        /*-------optimization problems where misconvergence occurs.-------------------------------*/
        if (strategy == 1) /* strategy DE0 (not in our paper) */
        {
            tmp = population[i];
            n = (int)(doubleRand(0,1)*D);
            L = 0;
            do
            {
                tmp[n] = bestit[n] + F*(population[r2][n]-population[r3][n]);
                n = (n+1)%D;
                L++;
            }while((doubleRand(0,1) < CR) && (L < D));
        }
            /*-------DE/rand/1/exp-------------------------------------------------------------------*/
            /*-------This is one of my favourite strategies. It works especially well when the-------*/
            /*-------"bestit[]"-schemes experience misconvergence. Try e.g. F=0.7 and CR=0.5---------*/
            /*-------as a first guess.---------------------------------------------------------------*/
        else if (strategy == 2) /* strategy DE1 in the techreport */
        {
            tmp = population[i];
            n = (int)(doubleRand(0,1)*D);
            L = 0;
            do
            {
                tmp[n] = population[r1][n] + F*(population[r2][n]-population[r3][n]);
                n = (n+1)%D;
                L++;
            }while((doubleRand(0,1) < CR) && (L < D));
        }
            /*-------DE/rand-to-best/1/exp-----------------------------------------------------------*/
            /*-------This strategy seems to be one of the best strategies. Try F=0.85 and CR=1.------*/
            /*-------If you get misconvergence try to increase NP. If this doesn't help you----------*/
            /*-------should play around with all three control variables.----------------------------*/
        else if (strategy == 3) /* similiar to DE2 but generally better */
        {
            tmp = population[i];
            n = (int)(doubleRand(0,1)*D);
            L = 0;
            do
            {
                tmp[n] = tmp[n] + F*(bestit[n] - tmp[n]) + F*(population[r1][n]-population[r2][n]);
                n = (n+1)%D;
                L++;
            }while((doubleRand(0,1) < CR) && (L < D));
        }
            /*-------DE/best/2/exp is another powerful strategy worth trying--------------------------*/
        else if (strategy == 4)
        {
            tmp = population[i];
            n = (int)(doubleRand(0,1)*D);
            L = 0;
            do
            {
                tmp[n] = bestit[n] +
                         (population[r1][n]+population[r2][n]-population[r3][n]-population[r4][n])*F;
                n = (n+1)%D;
                L++;
            }while((doubleRand(0,1) < CR) && (L < D));
        }
            /*-------DE/rand/2/exp seems to be a robust optimizer for many functions-------------------*/
        else if (strategy == 5)
        {
            tmp = population[i];
            n = (int)(doubleRand(0,1)*D);
            L = 0;
            do
            {
                tmp[n] = population[r5][n] +
                         (population[r1][n]+population[r2][n]-population[r3][n]-population[r4][n])*F;
                n = (n+1)%D;
                L++;
            }while((doubleRand(0,1) < CR) && (L < D));
        }

            /*=======Essentially same strategies but BINOMIAL CROSSOVER===============================*/

            /*-------DE/best/1/bin--------------------------------------------------------------------*/
        else if (strategy == 6)
        {
            tmp = population[i];
            n = (int)(doubleRand(0,1)*D);
            for (L=0; L<D; L++) /* perform D binomial trials */
            {
                if ((doubleRand(0,1) < CR) || L == (D-1)) /* change at least one parameter */
                {
                    tmp[n] = bestit[n] + F*(population[r2][n]-population[r3][n]);
                }
                n = (n+1)%D;
            }
        }
            /*-------DE/rand/1/bin-------------------------------------------------------------------*/
        else if (strategy == 7)
        {
            tmp = population[i];
            n = (int)(doubleRand(0,1)*D);
            for (L=0; L<D; L++) /* perform D binomial trials */
            {
                if ((doubleRand(0,1) < CR) || L == (D-1)) /* change at least one parameter */
                {
                    tmp[n] = population[r1][n] + F*(population[r2][n]-population[r3][n]);
                }
                n = (n+1)%D;
            }
        }
            /*-------DE/rand-to-best/1/bin-----------------------------------------------------------*/
        else if (strategy == 8)
        {
            tmp = population[i];
            n = (int)(doubleRand(0,1)*D);
            for (L=0; L<D; L++) /* perform D binomial trials */
            {
                if ((doubleRand(0,1) < CR) || L == (D-1)) /* change at least one parameter */
                {
                    tmp[n] = tmp[n] + F*(bestit[n] - tmp[n]) + F*(population[r1][n]-population[r2][n]);
                }
                n = (n+1)%D;
            }
        }
            /*-------DE/best/2/bin--------------------------------------------------------------------*/
        else if (strategy == 9)
        {
            tmp = population[i];
            n = (int)(doubleRand(0,1)*D);
            for (L=0; L<D; L++) /* perform D binomial trials */
            {
                if ((doubleRand(0,1) < CR) || L == (D-1)) /* change at least one parameter */
                {
                    tmp[n] = bestit[n] +
                             (population[r1][n]+population[r2][n]-population[r3][n]-population[r4][n])*F;
                }
                n = (n+1)%D;
            }
        }
            /*-------DE/rand/2/bin--------------------------------------------------------------------*/
        else
        {
            tmp = population[i];
            n = (int)(doubleRand(0,1)*D);
            for (L=0; L<D; L++) /* perform D binomial trials */
            {
                if ((doubleRand(0,1) < CR) || L == (D-1)) /* change at least one parameter */
                {
                    tmp[n] = population[r5][n] +
                             (population[r1][n]+population[r2][n]-population[r3][n]-population[r4][n])*F;
                }
                n = (n+1)%D;
            }
        }
    }

}
