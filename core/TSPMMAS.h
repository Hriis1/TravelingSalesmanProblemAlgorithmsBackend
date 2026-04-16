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
    struct Ant {
        std::vector<int> path;
        std::vector<bool> visited; 
        int currentCityIdx;  

        void reset(int N, int startCity) {
            path.clear();
            path.reserve(N);

            visited.assign(N, false);

            currentCityIdx = startCity;
            path.push_back(startCity);
            visited[startCity] = true;
        }
    };

public:

private:
};
