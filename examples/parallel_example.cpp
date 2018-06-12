/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde and Authors------------
----------------- e-mail: amin.kamal.2016@uni.strath.ac.uk ---------------------
--------------------------- Author: Amin Kamal ---------------------------------
*/

#include "stdio.h"
#include "Parallel/EvolutionaryAlgorithms.h"
#include "Parallel/PoolEA.h"
#include "Parallel/IMEA.h"
#include "Parallel/PSO.h"
#include "Parallel/GA.h"

using namespace Parallel;
using namespace std::chrono;

void PoolEATest()
{
    std::cout << "Pool-based:" << std::endl << std::endl;

    ///////////////////////////////////////////////////////////////////////////
    //Options
    //Problem options
    int n = 50;       //Problem dimension
    std::vector<std::vector<double>> bounds = {{-100, 100}};
    int func_num = 1; //Enter 0 to use the basic 1-d or 2-d problems defined in evaluateNonCEC() below,
    // enter > 0 for CEC2017 problems

    //Pool-based EA options
    int numberOfWorkers = 3;
    int populationSize = 100;
    int migrationRate = 250;
    int sampleSize = 30;
    int totalNumberOfGenerations = (n * 10000) / populationSize;
    int maxGenerationsPerWorker = totalNumberOfGenerations / numberOfWorkers;
    int initialPoolPopulationSize = populationSize * numberOfWorkers;
    auto em = ExecutionMethod::Async;

    //Output options
    int numberOfRuns = 3;
    bool printOutput = true;

    //Initialise EvoSpaceCPP container
    EvoSpaceCPP espp(initialPoolPopulationSize, n, bounds,
                     sampleSize, totalNumberOfGenerations, numberOfRuns, migrationRate, func_num);

    //Define algorithms
    
    // DE
    int crossOverStrategy = 3;
    double crossOverRatioDE = 1;
    double F_DE = 0.85;
    DifferentialEvolution de1(n, maxGenerationsPerWorker, crossOverRatioDE, crossOverStrategy, F_DE, populationSize, bounds, func_num);
    DifferentialEvolution de2(n, maxGenerationsPerWorker, crossOverRatioDE, crossOverStrategy, F_DE, populationSize, bounds, func_num);
    DifferentialEvolution de3(n, maxGenerationsPerWorker, crossOverRatioDE, crossOverStrategy, F_DE, populationSize, bounds, func_num);

    // PSO
    int numberOfAgents = 100;
    double F_PSO = 0.8;
    PSO pso1(n, populationSize, numberOfAgents, bounds, maxGenerationsPerWorker, F_PSO, func_num);
    PSO pso2(n, populationSize, numberOfAgents, bounds, maxGenerationsPerWorker, F_PSO, func_num);
    PSO pso3(n, populationSize, numberOfAgents, bounds, maxGenerationsPerWorker, F_PSO, func_num);

    // GA
    double crossOverRatioGA = 0.8;
    double mutationRatio = 0.15;
    GA ga1(n, populationSize, crossOverRatioGA, mutationRatio, maxGenerationsPerWorker, bounds, func_num);
    GA ga2(n, populationSize, crossOverRatioGA, mutationRatio, maxGenerationsPerWorker, bounds, func_num);
    GA ga3(n, populationSize, crossOverRatioGA, mutationRatio, maxGenerationsPerWorker, bounds, func_num);

    //Change algorithms and the number of workers here
    EvoSpaceWorker esw1(de1, espp);
    EvoSpaceWorker esw2(de2, espp);
    EvoSpaceWorker esw3(de3, espp);

    //Don't forget to add extra workers to the container (if any)
    espp.addWorker(esw1);
    espp.addWorker(esw2);
    espp.addWorker(esw3);

    ///////////////////////////////////////////////////////////////////////////
    //Output
    std::string em_string;
    if(em == ExecutionMethod::Sync) em_string = "Sync"; else em_string = "Async";
    //Generate filename
    std::string filename = "";
    if(func_num != 0)
        filename += "PoolEA_" + em_string + "_" + std::to_string(numberOfWorkers) + "_Workers_Problem_"
                           + std::to_string(func_num) + "_D" + std::to_string(n);
    else
        filename += "User-defined_Output";

    std::string constrains;
    bool constrains_bool = true;
    for(auto &b : bounds){
        constrains += "[";
        for(auto &bb : b){
            if(constrains_bool){
                constrains += std::to_string(bb) + ", ";
                constrains_bool = false;
            }
            else
                constrains += std::to_string(bb);
        }
        constrains += "]";
        constrains_bool = true;
    }
    std::string problem_type = "CEC Problem No." + std::to_string(func_num);
    if(func_num == 0)
        problem_type = "User-defined problem";
    char formatted_header[5000];
    const char * header =
            "\nStarting pool-based example with %d workers.\n"
            "Each worker will be filled with %d individuals.\n"
            "Each worker will then do %d generations\n"
            "while putting and taking %d best individuals\n"
            "from the pool every %d generations.\n\n"
            "Parameters:\n"
            "Problem type: %s\n"
            "Problem dimension: %d\n"
            "Problem constrains: %s\n"
            "Number of runs: %d\n"
            "Output filename: %s\n"
            "Print to stdout? %s\n";
    sprintf(formatted_header, header, numberOfWorkers, populationSize, maxGenerationsPerWorker, sampleSize,
            migrationRate, problem_type.c_str(), n, constrains.c_str(), numberOfRuns, filename.c_str(),
            printOutput ? "true" : "false");
    std::cout << formatted_header;
    //Write the above text to output file
    std::ofstream outfile;
    outfile.open ("output/" + filename + ".txt", std::ios_base::app);
    if(!outfile.good()) {
        std::cout << "OUTPUT ERROR: Ensure that a folder named 'output' exists!" << std::endl;
        return;
    }
    outfile << formatted_header;
    outfile.close();
    //Start the process
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    espp.start(em, filename, printOutput);
    //Store the time
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count();
    std::cout << std::endl << "The process took " << duration / 1e6 << "s to finish" << std::endl;
    outfile.open ("output/" + filename + ".txt", std::ios_base::app);
    if(!outfile.good()) {
        std::cout << "OUTPUT ERROR: Ensure that a folder named 'output' exists!" << std::endl;
        return;
    }
    outfile << std::endl << "The process took " << duration / 1e6 << "s to finish" << std::endl;
    outfile.close();
}

void IMEATest()
{
    std::cout << "Island-based:" << std::endl << std::endl;

    ///////////////////////////////////////////////////////////////////////////
    //Options

    //Problem options
    int n = 50;       //Problem dimension
    std::vector<std::vector<double>> bounds = {{-100, 100}};
    int func_num = 1; //Enter 0 to use the basic 1-d or 2-d problems defined in evaluateNonCEC() below,
    // enter > 0 for CEC2017 problems

    //Island-based EA options
    int numberOfIslands = 4;
    int demeSize = 30;
    int populationSize = 100;
    int migrationRate = 250;
    int totalNumberOfGenerations = (n * 10000) / populationSize;
    int maxGenerationsPerWorker = totalNumberOfGenerations / numberOfIslands;
    int initialPoolPopulationSize = populationSize * numberOfIslands;
    auto em = ExecutionMethod::Async;
    auto topology = Topology::Ring;

    //Output options
    int numberOfRuns = 3;
    bool printOutput = true;

    //Initialise Islands container
    IMEA islandmodelEA(demeSize, totalNumberOfGenerations, numberOfRuns, migrationRate, func_num);
    islandmodelEA.setTopology(topology, numberOfIslands);

    //Define algorithms

    // DE
    int crossOverStrategy = 3;
    double crossOverRatioDE = 1;
    double F_DE = 0.85;

    DifferentialEvolution de1(n, maxGenerationsPerWorker, crossOverRatioDE, crossOverStrategy, F_DE, populationSize, bounds, func_num);
    DifferentialEvolution de2(n, maxGenerationsPerWorker, crossOverRatioDE, crossOverStrategy, F_DE, populationSize, bounds, func_num);
    DifferentialEvolution de3(n, maxGenerationsPerWorker, crossOverRatioDE, crossOverStrategy, F_DE, populationSize, bounds, func_num);
    DifferentialEvolution de4(n, maxGenerationsPerWorker, crossOverRatioDE, crossOverStrategy, F_DE, populationSize, bounds, func_num);
  
    // PSO
    int numberOfAgents = 100;
    double F_PSO = 0.8;
    PSO pso1(n, populationSize, numberOfAgents, bounds, maxGenerationsPerWorker, F_PSO, func_num);
    PSO pso2(n, populationSize, numberOfAgents, bounds, maxGenerationsPerWorker, F_PSO, func_num);
    PSO pso3(n, populationSize, numberOfAgents, bounds, maxGenerationsPerWorker, F_PSO, func_num);
    PSO pso4(n, populationSize, numberOfAgents, bounds, maxGenerationsPerWorker, F_PSO, func_num);
  
    // GA
    double crossOverRatioGA = 0.8;
    double mutationRatio = 0.15;
    GA ga1(n, populationSize, crossOverRatioGA, mutationRatio, maxGenerationsPerWorker, bounds, func_num);
    GA ga2(n, populationSize, crossOverRatioGA, mutationRatio, maxGenerationsPerWorker, bounds, func_num);
    GA ga3(n, populationSize, crossOverRatioGA, mutationRatio, maxGenerationsPerWorker, bounds, func_num);
    GA ga4(n, populationSize, crossOverRatioGA, mutationRatio, maxGenerationsPerWorker, bounds, func_num);

    //Change algorithms and number of workers (islands) here
    Island isl1(1, islandmodelEA, de1);
    Island isl2(2, islandmodelEA, de2);
    Island isl3(3, islandmodelEA, de3);
    Island isl4(4, islandmodelEA, de4);

    //Don't forget to add extra islands to the container (if any)
    islandmodelEA.addIsland(isl1);
    islandmodelEA.addIsland(isl2);
    islandmodelEA.addIsland(isl3);
    islandmodelEA.addIsland(isl4);


    ///////////////////////////////////////////////////////////////////////////
    //Output
    std::string topology_string;
    if(topology == Topology::Ring) topology_string = "Ring";
    else if(topology == Topology::FullyConnected) topology_string = "FullyConnected";
    else if(topology == Topology::Broadcast) topology_string = "Broadcast";
    std::string em_string;
    if(em == ExecutionMethod::Sync) em_string = "Sync"; else em_string = "Async";
    //Generate filename
    std::string filename = "";
    if(func_num != 0)
        filename += "IMEA_" + topology_string + "_" + em_string + "_" + std::to_string(numberOfIslands) + "_Islands_Problem_"
                    + std::to_string(func_num) + "_D" + std::to_string(n);
    else
        filename += "User-defined_Output";

    std::string constrains;
    bool constrains_bool = true;
    for(auto &b : bounds){
        constrains += "[";
        for(auto &bb : b){
            if(constrains_bool){
                constrains += std::to_string(bb) + ", ";
                constrains_bool = false;
            }
            else
                constrains += std::to_string(bb);
        }
        constrains += "]";
        constrains_bool = true;
    }
    std::string problem_type = "CEC Problem No." + std::to_string(func_num);
    if(func_num == 0)
        problem_type = "User-defined problem";
    char formatted_header[5000];
    const char * header =
            "\nStarting island-based example with %d islands.\n"
                    "\nSelected topology is %s.\n"
                    "Each island will start with %d individuals.\n"
                    "Each island will then do %d generations\n"
                    "while sending and accepting %d best individuals\n"
                    "to and from its neighbours every %d generations.\n\n"
                    "Parameters:\n"
                    "Problem type: %s\n"
                    "Problem dimension: %d\n"
                    "Problem constrains: %s\n"
                    "Number of runs: %d\n"
                    "Output filename: %s\n"
                    "Print to stdout? %s\n";
    sprintf(formatted_header, header, numberOfIslands, topology_string.c_str(), populationSize, maxGenerationsPerWorker, demeSize,
            migrationRate, problem_type.c_str(), n, constrains.c_str(), numberOfRuns, filename.c_str(),
            printOutput ? "true" : "false");
    std::cout << formatted_header;
    //Write the above text to output file
    std::ofstream outfile;
    outfile.open ("output/" + filename + ".txt", std::ios_base::app);
    if(!outfile.good()) {
        std::cout << "OUTPUT ERROR: Ensure that a folder named 'output' exists!" << std::endl;
        return;
    }
    outfile << formatted_header;
    outfile.close();
    //Start the process
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    islandmodelEA.start(em, filename, printOutput);
    //Store the time
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count();
    std::cout << std::endl << "The process took " << duration / 1e6 << "s to finish" << std::endl;
    outfile.open ("output/" + filename + ".txt", std::ios_base::app);
    if(!outfile.good()) {
        std::cout << "OUTPUT ERROR: Ensure that a folder named 'output' exists!" << std::endl;
        return;
    }
    outfile << std::endl << "The process took " << duration / 1e6 << "s to finish" << std::endl;
    outfile.close();
}

void Parallel::evaluateNonCEC(Individual &ind)
{
    double x1,x2;
    x1 = ind.getMembers()[0]; //0 since 1-D
    x2 = ind.getMembers()[1]; //1 since 2-D
    //ind.setFitness(x1*x1-4*x1+2); //1-D problem
    ind.setFitness(x2*x2+x1*x2+x1*x1-4*x1+2); //2-D problem
}

int main(){
    std::cout << "Welcome to SMART-O2CPP2 Parallel Example!" << std::endl << std::endl;

    PoolEATest(); // Remember to select the required test!
    //IMEATest();

    return 0;
}

