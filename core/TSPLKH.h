#pragma once
#include <climits>
#include <array>
#include <cassert>

#include "TSPAlgo.h"
#include "TSPSolution.h"
#include "TSPUtils.h"


struct LKHConfig
{
    int maxTrials = 0;
    int maxCandidates = 0;
    int maxDepth = 0;
    int backtrackingLimit = 0;
    int runs = 0;
    int kickStrength = 0;

    bool useAlphaCandidates = false;
    bool useDontLookBits = false;
    bool useBacktracking = false;
};


class TSPMMAS: public TSPAlgo
{
private:

public:

    TSPMMAS(LKHConfig& config, unsigned int seed = std::random_device{}())
        :TSPAlgo(seed), _config(config)
    {}

    void solve(const std::vector<std::vector<int>>& adjMat) override
    {
    }

private:

    LKHConfig _config;
};
