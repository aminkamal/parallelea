/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
------------Copyright (C) 2017 University of Strathclyde and Authors------------
----------------- e-mail: amin.kamal.2016@uni.strath.ac.uk ---------------------
--------------------------- Author: Amin Kamal ---------------------------------
*/

#ifndef POPULATION_H
#define POPULATION_H

#include <vector>
#include <deque>
#include <algorithm>

namespace Parallel
{
    class Individual
    {
    private:
        int id_;
        std::vector<double> members_;
        double fitness_;

    public:
        //get/set
        int getID() const { return id_; };
        void setID(int id) { id_ = id; };
        std::vector<double> getMembers() const { return members_; };
        void setMembers(std::vector<double> fitness) { members_ = fitness; };
        double getFitness() const { return fitness_; };
        void setFitness(double fitness) { fitness_ = fitness; };

        Individual(int id, std::vector<double> members, double fitness)
                : id_(id), members_(members), fitness_(fitness) {}

        double& operator[] (int i) {
            return members_[i];
        }

        void operator= (std::vector<double> i) {
            members_ = i;
        }

        Individual& operator+(Individual lhs)
        {
            for(int j=0;j<members_.size();j++)
                members_[j] += lhs[j];
            return *this;
        }

        bool operator<(double const& b) const
        {
            return fitness_ < b;
        }

        bool operator<(Individual const& b) const
        {
            return fitness_ < b.getFitness();
        }

        std::string toString(int func_num = 0)
        {
            std::string temp = "Individual ID: " + std::to_string(id_);

            if(func_num == 0){
                temp += ", Fitness = ";
                temp += std::to_string(fitness_);
            }
            else{
                temp += ", Error = ";
                temp += std::to_string( std::abs(fitness_ - (func_num * 100 )) );
            }
            temp += ", Members {";
            for(int j=0;j<members_.size();j++)
                temp += std::to_string(members_[j]) + ", ";
            temp += "}";
            return temp;
        }
    };

    class Population : public std::deque<Individual>
    {
    private:
        int id_;
    public:
        int getID() const { return id_; };
        void setID(int id) { id_ = id; };
        void sort()
        {
            std::sort(begin(), end()); //the underlying inherited deque iterator
        }
    };
}

#endif //SMART_O2CPP_POPULATION_H
