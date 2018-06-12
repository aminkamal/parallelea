/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde and Authors------------
----------------- e-mail: amin.kamal.2016@uni.strath.ac.uk ---------------------
--------------------------- Author: Amin Kamal ---------------------------------
*/

#include "Parallel/IMEA.h"

namespace Parallel
{
    //
    // Island code
    //

    bool Island::isFinished()
    {
        if(ea_.getNumberOfGenerations() < imea_.getMaxNumberOfGenerations())
            return false;
        ea_.incrementRunCount();
        return true;
    }

    int Island::getID()
    {
        return islandID_;
    }

    void Island::initEA()
    {
        buffer.clear();
        ea_.generatePopulationForIMEA();
        ea_.init();
    }

    void Island::appendToBuffer(Population p)
    {
        buffer.insert(buffer.end(), p.begin(), p.end());
    }

    void Island::processBuffer()
    {
        if(!buffer.empty())
        {
            ea_.acceptBuffer(buffer);
            buffer.clear();
        }
    }

    void Island::executeSync()
    {
        int i = 0;
        while(!isFinished())
        {
            if (i > imea_.getMigrationRate()) {
                return;
            }else{
                ea_.executeEA();
            }
            i++;
        }
        imea_.syncFinishCounter++;
    }

    void Island::executeAsync()
    {
        int i = 0;
        while(!isFinished())
        {
            if (i > imea_.getMigrationRate()) {
                imea_.islandMutex.lock();
                //Migrate to all neighbours
                imea_.spreadPopulation(
                        getID(),
                        ea_.getIndividuals(imea_.getDemeSize(), false)
                );
                processBuffer();
                i = 0;
                imea_.islandMutex.unlock();
            }else{
                ea_.executeEA();
            }
            i++;
        }
    }

    //
    // IMEA code
    //

    void IMEA::start(ExecutionMethod em, std::string filename, bool printError)
    {
#ifndef USEMULTITHREADS
        std::cout << "Warning: The program was not compiled with threading support." << std::endl;
        std::cout << "Warning: Ignoring Sync and Async execution methods directives." << std::endl;
#endif
        std::thread threads[maxNumberOfThreads];
        int numberOfIslands = (int) islV.size();

        for(int i = 0; i<numberOfIslands;i++)
            islV[i].initEA();

        for(int k=0;k<numberOfRuns_;k++){
#ifdef USEMULTITHREADS
            if(em == ExecutionMethod::Sync){
                syncFinishCounter = 0;
                while(syncFinishCounter < islV.size())
                {
                    for (int i=0; i<numberOfIslands; i++){
                        threads[i] = std::thread(&Island::executeSync, islV[i]);
                    }
                    for (int i=0; i<numberOfIslands; i++){
                        threads[i].join();
                    }
                    for (int i=0; i<numberOfIslands; i++){
                        //Migrate to all neighbours
                        spreadPopulation(
                                islV[i].getID(),
                                islV[i].ea_.getIndividuals(demeSize_, false)
                        );
                    }
                    for (int i=0; i<numberOfIslands; i++){
                        islV[i].processBuffer();
                    }
                }
            }else if(em == ExecutionMethod::Async){

                for (int i=0; i<numberOfIslands; i++){
                    threads[i] = std::thread(&Island::executeAsync, islV[i]);
                }
                for (int i=0; i<numberOfIslands; i++){
                    threads[i].join();
                }
            }
#else
            std::vector<int> j;
            j.reserve(sizeof(islV));
            int currentlyProcessing = 0, finishedCount = 0;
            while(finishedCount < numberOfIslands){
                islV[currentlyProcessing].ea_.executeEA();

                if(islV[currentlyProcessing].ea_.getNumberOfGenerations() == migrationRate_){
                    //Migrate to all neighbours
                    spreadPopulation(
                            islV[currentlyProcessing].getID(),
                            islV[currentlyProcessing].ea_.getIndividuals(demeSize_, false)
                    );
                    islV[currentlyProcessing].processBuffer();
                    j[currentlyProcessing] = 0;
                }

                if(islV[currentlyProcessing].isFinished())
                    finishedCount++;

                j[currentlyProcessing]++;
                currentlyProcessing++;

                if(currentlyProcessing == numberOfIslands)
                    currentlyProcessing = 0;
            }
#endif

            if(k+1 == numberOfRuns_)
            for(int i=0; i<numberOfIslands;i++)
            {
                if(func_num_ != 0){
                    if(printError)
                        islV[i].ea_.printCECError();
                    islV[i].ea_.saveOutputToFile("Island" + std::to_string(i + 1), filename);
                    islV[i].ea_.saveCECErrorToFile("Island" + std::to_string(i + 1), filename);
                }else{
                    if(printError)
                        islV[i].ea_.output();
                    islV[i].ea_.saveOutputToFile("Island" + std::to_string(i + 1), filename);
                }
            }

            //Reset
            if(func_num_ != 0)
            for(int i=0; i<numberOfIslands;i++)
            {
                islV[i].ea_.reset();
                islV[i].initEA();
            }
        }
    }

    void IMEA::addIsland(Island isl)
    {
        islV.push_back(isl);
    }

    int IMEA::getMaxNumberOfGenerations()
    {
        return workerGenerations_ / (int)islV.size();
    }

    int IMEA::getMigrationRate()
    {
        return migrationRate_;
    }

    int IMEA::getDemeSize()
    {
        return demeSize_;
    }

    void IMEA::spreadPopulation(int islandID, Population p)
    {
        for(auto i : g.getNeighbours(islandID)){
            islV[i-1].appendToBuffer(p); //Our code starts with index 1
        }
    }

    void IMEA::setTopology(Topology t, int numberOfIslands)
    {
        switch(t)
        {
            case Topology::Ring:
                /*
                 RING TOPOLOGY
                 Example of generated graph if numberOfIslands = 8
                 1 -> {2, 8}
                 2 -> {1, 3}
                 3 -> {2, 4}
                 4 -> {3, 5}
                 5 -> {4, 6}
                 6 -> {5, 7}
                 7 -> {6, 8}
                 8 -> {7, 1}
                 */
                //Connect island number 1 with the last island
                g.addNode(1);
                g.addNeighbour(1, numberOfIslands);
                g.addNeighbour(1, 2);
                for(int i=2;i<=numberOfIslands;i++){
                    g.addNode(i);
                    g.addNeighbour(i, i-1);
                    if(i != numberOfIslands)
                        g.addNeighbour(i, i+1);
                }
                g.addNeighbour(numberOfIslands, 1);
                break;
            case Topology::FullyConnected:
                 /*
                 FULLY CONNECTED TOPOLOGY
                 1 -> {2, 3, 4, 5, 6, 7, 8}
                 2 -> {1, 3, 4, 5, 6, 7 ,8}
                 3 -> {1, 2, 4, 5, 6, 7, 8}
                 4 -> {1, 2, 3, 5, 6, 7, 8}
                 5 -> {1, 2, 3, 4, 6, 7, 8}
                 6 -> {1, 2, 3, 4, 5, 7, 8}
                 7 -> {1, 2, 3, 4, 5, 6, 8}
                 8 -> {1, 2, 3, 4, 5, 6, 7}
                 */
                for(int i=1;i<=numberOfIslands;i++){
                    g.addNode(i);
                    for(int j=1;j<=numberOfIslands;j++)
                        if(i != j)
                            g.addNeighbour(i, j);
                }
                break;
            case Topology::Broadcast:
                /*
                 BROADCAST TOPOLOGY
                 1 -> {2, 3, 4, 5, 6, 7, 8};
                 2 -> {1};
                 3 -> {1};
                 4 -> {1};
                 5 -> {1};
                 6 -> {1};
                 7 -> {1};
                 8 -> {1};
                */
                g.addNode(1);
                for(int i=2;i<=numberOfIslands;i++) {
                    g.addNeighbour(1, i);

                    g.addNode(i);
                    g.addNeighbour(i, 1);
                }
                break;
        }
    }
}
