#pragma once
#include <climits>
#include <array>
#include <cassert>

#include "TSPAlgo.h"
#include "TSPSolution.h"
#include "TSPUtils.h"


struct LKHConfig
{
    // ASIGNED VALUES ARE JUST DEFAULTS

    int maxTrials = 50;                 //How many attempts LKH makes
    int maxCandidates = 15;             //How many candidate edges each city considers
    int maxDepth = 3;                   //Maximum depth of the variable k-opt search
    int backtrackingLimit = 20;         //Limits how many failed alternatives are explored.
    int runs = 1;                       //How many independent full runs to do
    int kickStrength = 4;               //How strong the perturbation is when stuck


    long long precision = 100;                //precision for calculating transformed costs
};


class TSPLKH: public TSPAlgo
{
private:

public:

    TSPLKH(const LKHConfig& config, unsigned int seed = std::random_device{}())
        :TSPAlgo(seed), _config(config)
    {}

    void solve(const std::vector<std::vector<int>>& adjMat) override
    {
        //num cities
        int n = adjMat.size();

        //asign penalties to 0 at the start
        pi.assign(n, 0);
    }

private:

    //Returns the transformed cost between 2 cities
    long long getTransformedCost(int i, int j, const std::vector<std::vector<int>>& adjMat) const
    {
        return (long long)(_config.precision * (long long)adjMat[i][j] + pi[i] + pi[j]);
    }

private:
    LKHConfig _config; //config data for solver

    std::vector<long long> pi; //pnalties

};
