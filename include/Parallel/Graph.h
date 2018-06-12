/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde and Authors------------
----------------- e-mail: amin.kamal.2016@uni.strath.ac.uk ---------------------
--------------------------- Author: Amin Kamal ---------------------------------
*/

#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <map>
#include <iostream>

namespace Parallel{
    class UndirectedGraph{
        private:
        std::map<int, std::vector<int>> graph;

        public:
        void addNode(int node_id){
            graph[node_id];
        }

        void addNeighbour(int node_id, int neighbour_id){
            graph[node_id].push_back(neighbour_id);
        }

        std::vector<int> getNeighbours(int node_id){
            return graph[node_id];
        }

        //Debug
        void printTopology(){
            for(unsigned int i=1;i<=graph.size();i++){
                std::cout << i << " -> {";
                for(unsigned int j=0;j<graph[i].size();j++)
                    std::cout << graph[i][j];
                std::cout << "} ";
            }
        }
    };
}

#endif