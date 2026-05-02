#pragma once
#include <climits>
#include <array>
#include <cassert>

#include "TSPAlgo.h"
#include "TSPSolution.h"
#include "TSPUtils.h"


class TSPMMAS: public TSPAlgo
{
private:

public:

    TSPMMAS(unsigned int seed = std::random_device{}())
        :TSPAlgo(seed)
    {}

    void solve(const std::vector<std::vector<int>>& adjMat) override
    {
    }

private:
};
