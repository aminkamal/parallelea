C++ Parallel island and pool-based evolutionary algorithms
==========================================================

This project was part of SMART-O2CPP (Strathclyde Mechanical and Aerospace Research Toolbox for Optimisation and Optimal Control), it was developed by Amin Kamal (me) as part of my MSc thesis titled: "On the comparison of parallel evolutionary algorithms" supervised by Dr. Annalisa Riccardi.

It was used to compare three evolutionary-based algorithms: Differential evolution (DE), Genetic Algorithms (GA) and Particle swarm optimisation (PSO) in two different parallelisation paradigms (Pool-based and Island-based).

View the commented example file /examples/Parallel.cpp to see how to use the library and to view supported parameters and options.

About `SMART`
------
`SMART` is a collection of toolboxes developed and maintained since 2015 by the department of Mechanical and Aerospace Engineer of Strathclyde University. `SMART-O2CPP` in particular, is a collection of optimisation and optimal control transcription techniques.

`SMART` is aimed at making code development inside the team more efficient, reusable and easy maintainabile.

`SMART` is CMake-based C++ project. It can be used as the basis for new projects. 

Features
------
* User-defined functions support; The project was originally intended to be used with the CEC2017 problem set, but simple 1D and 2D and n-D single objective functions such as: y=x1*x1-4*x1+2 and y=x2*x2+x1*x2+x1*x1-4*x1+2 are also supported

* Extensibility; although it was initially designed with 3 algorithms in mind, it is easy to add more algorithms thanks to the extensible code architecture

* Easy to read; you don't need a PhD in Computer Science and 10 years of C++ voodoo experience to understand and modify the source code

Requirements
------
* C++11 compatible compiler
  GCC prefered
  Intel ICC works
  MSVC might work

* CMake 3.0+

Installation
------

Run the following commands to download, build, and install this project.
* git clone this repo
* cd smart-o2cpp
* edit the file examples/Parallel.cpp with your algorithms, settings..etc
* mkdir build
* mkdir build/output
* cp -R cec2017_bound_constrained_input_data build/cec2017_bound_constrained_input_data
* cd build && cmake .. && make
* ./parallel

License
------

Mozilla Public License Version 2.0

Disclaimer
------

The copyright holders are not liable for any damage(s) incurred due to improper use of `smart-o2cpp`.

Note on automation of runs
------

The example supplied is a simplified one. It is possible to automate the runs instead of manually editing and recompiling the project everytime, and this is what was done during my project. For example, the following code should get you started on how to automate runs:
```C++
    int n[4] = {10, 30, 50, 100};
    int n_mig[4] = {50, 150, 250, 500};
    std::string combinations[15] = {
            "0000",
            "0001",
            "0002",
            "0011",
            "0012",
            "0022",
            "0111",
			........
    };
    EvolutionaryAlgorithm** eas = new EvolutionaryAlgorithm*[4];

    for(int ddd = 0;ddd<4;ddd++)
        for(int method = 0; method < 2;method++){
            if(method == 0) em = ExecutionMethod::Async;
            if(method == 1) em = ExecutionMethod::Sync;

            int totalNumberOfGenerations = (n[ddd] * 10000) / populationSize;
            int maxGenerationsPerWorker = totalNumberOfGenerations / numberOfWorkers;

            for(int func_num=1;func_num<31;func_num++)
                for(auto &c : combinations){
                    int deidx=0;
                    int psoidx=0;
                    int gaidx=0;

                    EvoSpaceCPP espp(initialPoolPopulationSize, n[ddd], bounds,
                                     sampleSize, totalNumberOfGenerations, numberOfRuns, n_mig[ddd], func_num);

                    DifferentialEvolution de1(n[ddd], maxGenerationsPerWorker, 1, 3, 0.85, populationSize, bounds, func_num);
                    DifferentialEvolution de2(n[ddd], maxGenerationsPerWorker, 1, 3, 0.85, populationSize, bounds, func_num);
					........

                    PSO pso1(n[ddd], populationSize, 100, bounds, maxGenerationsPerWorker, 0.8 , func_num);
                    PSO pso2(n[ddd], populationSize, 100, bounds, maxGenerationsPerWorker, 0.8 , func_num);
					........

                    GA ga1(n[ddd], populationSize, 0.8, 0.15, maxGenerationsPerWorker, bounds, func_num);
                    GA ga2(n[ddd], populationSize, 0.8, 0.15, maxGenerationsPerWorker, bounds, func_num);
					........

                    for(int i=0;i<numberOfWorkers;i++){
                        if(c[i] == '0'){
                            if(deidx == 0) eas[i] = &de1;
                            else if(deidx == 1) eas[i] = &de2;
                            else if(deidx == 2) eas[i] = &de3;
                            else if(deidx == 3) eas[i] = &de4;
                            deidx++;
                        }else if(c[i] == '1'){
                            if(psoidx == 0) eas[i] = &pso1;
                            else if(psoidx == 1) eas[i] = &pso2;
                            else if(psoidx == 2) eas[i] = &pso3;
                            else if(psoidx == 3) eas[i] = &pso4;
                            psoidx++;
                        }else if(c[i] == '2'){
                            if(gaidx == 0) eas[i] = &ga1;
                            else if(gaidx == 1) eas[i] = &ga2;
                            else if(gaidx == 2) eas[i] = &ga3;
                            else if(gaidx == 3) eas[i] = &ga4;
                            gaidx++;
                        }
                    }

                    EvoSpaceWorker esw1(*eas[0], espp);
                    EvoSpaceWorker esw2(*eas[1], espp);
					........

                    espp.addWorker(esw1);
                    espp.addWorker(esw2);
					........
```

(The code used in the thesis which automated 10800 runs is available on request).
