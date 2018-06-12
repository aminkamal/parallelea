/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde and Authors------------
----------------- e-mail: amin.kamal.2016@uni.strath.ac.uk ---------------------
--------------------------- Author: Amin Kamal ---------------------------------
*/

#include "Parallel/PoolEA.h"

namespace Parallel
{

    //
    // EvoSpaceWorker code
    //

    void EvoSpaceWorker::take(int n)
    {
        ea_.acceptTake(
                espp_.get_sample(n)
        );
    }

    void EvoSpaceWorker::storePopulationIntoPool(int n)
    {
        espp_.put_sample(
                ea_.getIndividuals(n)
        );
    }

    bool EvoSpaceWorker::isFinished()
    {
        if(ea_.getNumberOfGenerations() < espp_.getMaxNumberOfGenerations())
            return false;
        ea_.incrementRunCount();
        return true;
    }

    void EvoSpaceWorker::executeSync()
    {
        int i = 0;
        while(!isFinished())
        {
            if (i == espp_.getMigrationRate()) {
                return;
            }else{
                ea_.executeEA();
            }
            i++;
        }
        espp_.syncFinishCounter++;
    }

    void EvoSpaceWorker::executeAsync()
    {
        int i = 0;
        while(!isFinished())
        {
            if (i == espp_.getMigrationRate()) {
                //std::vector is not thread safe
                espp_.poolMutex.lock();
                storePopulationIntoPool(espp_.getSampleSize());
                take(espp_.getSampleSize());
                i = 0;
                espp_.poolMutex.unlock();
            }else{
                ea_.executeEA();
            }
            i++;
        }
    }

    //
    // EvoSpaceCPP (the container/pool) code
    //

    void EvoSpaceCPP::start(ExecutionMethod em, std::string filename, bool printError)
    {
#ifndef USEMULTITHREADS
        std::cout << "Warning: The program was not compiled with threading support." << std::endl;
        std::cout << "Warning: Ignoring Sync and Async execution methods directives." << std::endl;
#endif
        std::thread threads[maxNumberOfThreads];
        int numberOfWorkers = (int) eswV.size();

        // Fill workers with initial population from the pool
        distributeFairly(populationSize_ / numberOfWorkers, eswV);

        for(int k=0;k<numberOfRuns_;k++){
#ifdef USEMULTITHREADS
            if(em == ExecutionMethod::Sync){
                syncFinishCounter = 0;
                while(syncFinishCounter < numberOfWorkers)
                {
                    for (int i=0; i<numberOfWorkers; i++){
                        threads[i] = std::thread(&EvoSpaceWorker::executeSync, eswV[i]);
                    }
                    for (int i=0; i<numberOfWorkers; i++){
                        threads[i].join();
                    }
                    for (int i=0; i<numberOfWorkers; i++){
                        eswV[i].storePopulationIntoPool(sampleSize_);
                    }
                    distributeFairly(sampleSize_, eswV);
                }
            }else if(em == ExecutionMethod::Async){
                for (int i=0; i<numberOfWorkers; i++) {
                    threads[i] = std::thread(&EvoSpaceWorker::executeAsync, eswV[i]);
                }
                for (int i=0; i<numberOfWorkers; i++){
                    threads[i].join();
                }
            }
#else
            std::vector<int> j;
            j.reserve(sizeof(eswV));
            int currentlyProcessing = 0, finishedCount = 0;
            while(finishedCount < numberOfWorkers){
                eswV[currentlyProcessing].ea_.executeEA();

                if(eswV[currentlyProcessing].ea_.getNumberOfGenerations() == migrationRate_){
                    eswV[currentlyProcessing].storePopulationIntoPool(sampleSize_);
                    eswV[currentlyProcessing].take(sampleSize_);
                    j[currentlyProcessing] = 0;
                }

                if(eswV[currentlyProcessing].isFinished())
                    finishedCount++;

                j[currentlyProcessing]++;
                currentlyProcessing++;

                if(currentlyProcessing == numberOfWorkers)
                    currentlyProcessing = 0;
            }
#endif

            // Combine all individuals from all workers into the main pool
            for(int i=0; i<numberOfWorkers;i++) {
                Population poptobecopied = eswV[i].ea_.getPopulation();
                pool_population.insert(pool_population.end(),
                                           poptobecopied.begin(),
                                           poptobecopied.end());
            }
            pool_population.sort();

            if(printError){
                std::cout << "========================" << std::endl;
                if(func_num_ == 0)
                    std::cout << "     Pool population    " << std::endl;
                else
                    std::cout << "Pool population (CEC Error values)" << std::endl;
                std::cout << "========================" << std::endl;
                for(auto &p : pool_population){
                    if(func_num_ == 0)
                        std::cout << p.toString() << std::endl;
                    else
                        std::cout << p.toString(func_num_) << std::endl;
                }
            }

            // Write pool population (the solution) to file
            std::ofstream outfile;
            outfile.open ("output/" + filename + ".txt", std::ios_base::app);
            if(!outfile.good()){
                std::cout << "OUTPUT ERROR: Ensure that a folder named 'output' exists!" << std::endl;
                return;
            }
            outfile << "========================" << std::endl;
            if(func_num_ == 0)
                outfile << "     Pool population    " << std::endl;
            else
                outfile << "Pool population (CEC Error values)" << std::endl;
            outfile << "========================" << std::endl;
            for(auto &p : pool_population){
                if(func_num_ == 0)
                    outfile << p.toString() << std::endl;
                else
                    outfile << p.toString(func_num_) << std::endl;
            }
            outfile << "\n";
            outfile.close();

            //is this the last iteration? if so, save results before resetting
            if(k+1 == numberOfRuns_){
                for(int i=0; i<numberOfWorkers;i++) {
                    if (func_num_ != 0) {
                        eswV[i].ea_.saveOutputToFile("Worker" + std::to_string(i + 1), filename);
                    } else {
                        if (printError)
                            eswV[i].ea_.output();
                    }
                }
            }

            //Reset
            if(func_num_ != 0)
            {
                pool_population.clear();
                sample_queue.clear();
                pool_population = pool_population = EvolutionaryAlgorithm::generatePopulation(D_, populationSize_,
                                                                                              bounds, func_num_);
                for(int i=0; i<numberOfWorkers;i++)
                    eswV[i].ea_.reset();
                distributeFairly(populationSize_ / numberOfWorkers, eswV);
            }
        }
    }

    void EvoSpaceCPP::addWorker(EvoSpaceWorker esw)
    {
        eswV.push_back(esw);
    }

    int EvoSpaceCPP::getMaxNumberOfGenerations()
    {
        return workerGenerations_ / (int)eswV.size();
    }

    int EvoSpaceCPP::getMigrationRate()
    {
        return migrationRate_;
    }

    int EvoSpaceCPP::getSampleSize()
    {
        return sampleSize_;
    }

    void EvoSpaceCPP::distributeFairly(int n, std::vector<EvoSpaceWorker>& eswV)
    {
        Population popToBeDistributed = get_sample(n*eswV.size());
        if(TAKET == TAKETYPE::BEST){
            // Distribute the best population from the pool randomly to each worker
            // to avoid giving the first worker advantage over others
            std::shuffle (popToBeDistributed.begin(),
                          popToBeDistributed.end(),
                          std::default_random_engine(EvolutionaryAlgorithm::s_seed));
        }
        for (int i=0; i<eswV.size(); i++){
            Population subPop;
            Population::iterator it = popToBeDistributed.begin();
            while (it != (popToBeDistributed.begin() + n))
                subPop.push_back(*it++);
            popToBeDistributed.erase(popToBeDistributed.begin(), popToBeDistributed.begin() + n);
            eswV[i].ea_.acceptTake(subPop);
        }
    }

    // The following code was converted from EvoSpacy.py

    Population EvoSpaceCPP::get_sample(int size)
    {
        if (pool_population.size() <= MIN_SIZE){
            respawn(RE_INSERT_SAMPLES);
        }
        sample_counter++;

        Population population;
        population.setID(sample_counter);

        //Get random individuals from the population (as per the specifications)
        if(TAKET == TAKETYPE::RANDOM){
            for(unsigned int i = 0;i<size;i++){
                int luckyWinner = (int)EvolutionaryAlgorithm::s_doubleRand(0, (int) pool_population.size() - 1);
                population.push_back(pool_population[luckyWinner]);
                pool_population.erase(pool_population.begin()+luckyWinner);
            }
        }else{
            pool_population.sort();
            for(unsigned int i = 0;i<size;i++){
                population.push_back(pool_population[i]);
            }
            pool_population.erase(pool_population.begin(), pool_population.begin() + size - 1);
        }

        sample_queue.push_back( population );
        return population;
    }

    void EvoSpaceCPP::put_sample(Population population)
    {
        returned_counter++;
        std::deque<Individual>::iterator it = population.begin();
        while (it != population.end()){
            pool_population.push_back(*it++);
        }

        std::deque<Population>::iterator it2 = sample_queue.begin();
        while (it2 != sample_queue.end()){
            if((*it2).getID() == population.getID()){
                sample_queue.erase(it2);
                break; //we only need to remove the first occurrence
            }
            it2++;
        }
    }

    void EvoSpaceCPP::respawn_sample(Population population)
    {
        std::vector<Individual> members;
        std::deque<Individual>::iterator it = population.begin();
        while (it != population.end()){
            pool_population.push_back(*it);
            it++;
        }
    }

    void EvoSpaceCPP::respawn(int n)
    {
        if (RESPAWN == Parallel::RESPAWNTYPE::REINSERT){
            int current_size = (int) sample_queue.size();
            int nn;
            if (n > current_size)
                nn = current_size;
            else
                nn = n;
            for(int i = 0;i<nn;i++){
                respawn_sample( sample_queue.front() );
                sample_queue.pop_front();
            }
        }else if(RESPAWN == Parallel::RESPAWNTYPE::RANDOM){
            Population randomPop = EvolutionaryAlgorithm::generatePopulation(D_, n, bounds, func_num_);
            pool_population.insert(pool_population.end(), randomPop.begin(), randomPop.end());
        }
    }
}
