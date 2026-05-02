#pragma once
#include <climits>
#include <array>
#include <cassert>

#include "TSPAlgo.h"
#include "TSPSolution.h"
#include "TSPUtils.h"


struct LKHConfig
{
    //Default values
    int maxTrials = 50;
    int maxCandidates = 15;
    int maxDepth = 3;
    int backtrackingLimit = 20;
    int runs = 1;
    int kickStrength = 4;
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
    }

private:

private:
    LKHConfig _config;
};
