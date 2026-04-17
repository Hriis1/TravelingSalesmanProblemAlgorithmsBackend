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

    TSPMMAS(const std::vector<std::vector<int>>& adjMat, int numAnts, int numIterations, double alpha, double beta, double rho, double tauMin, double tauMax, unsigned int seed = std::random_device{}()) 
        :TSPAlgo(adjMat, seed), _numAnts(numAnts), _numIterations(numIterations), _alpha(alpha), _beta(beta), _rho(rho), _tauMin(tauMin), _tauMax(tauMax)
    {}

private:
    int _numAnts;        // number of ants per iteration
    int _numIterations;  // total iterations

    double _alpha;       // pheromone importance
    double _beta;        // heuristic importance
    double _rho;         // evaporation rate (0–1)

    //Pheromone bounds
    double _tauMin;
    double _tauMax;

    int _globalBestFreq = 10;

    std::vector<Ant> _ants;
    std::vector<std::vector<double>> _pheromone;
    std::vector<std::vector<double>> _heuristic;
};
