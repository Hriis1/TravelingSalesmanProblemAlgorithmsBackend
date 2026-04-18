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

    TSPMMAS(int numAnts, int numIterations, float alpha, float beta, float rho, unsigned int seed = std::random_device{}()) 
        :TSPAlgo(seed), _numAnts(numAnts), _ants(numAnts), _numIterations(numIterations), _alpha(alpha), _beta(beta), _rho(rho)
    {}

    void solve(const std::vector<std::vector<int>>& adjMat) override
    {
        int nCities = adjMat.size();

        //Initial solution using nearst neighbor
        _currSolution.path = TSPUtils::nearestNeighborPath(adjMat);
        _currSolution.calculateDist(adjMat);

        //Init _tauMin and _tauMax
        _tauMax = 1.0f / (_rho * _currSolution.dist);

        float pBest = 0.05f;
        float pRoot = pow(pBest, 1.0f / nCities);
        float avg = nCities / 2.0f;

        _tauMin = _tauMax * (1 - pRoot) / ((avg - 1) * pRoot);

        //Init pheromone matrix with _tauMax
        _pheromone.assign(nCities, std::vector<float>(nCities, _tauMax));

        //Init the heuristic matrix with 1/<dist between cities>
        _heuristic.assign(nCities, std::vector<float>(nCities));
        for (size_t y = 0; y < nCities; y++)
        {
            for (size_t x = 0; x < nCities; x++)
            {
                if (y != x && adjMat[y][x] > 0)
                    _heuristic[y][x] = 1.0f / adjMat[y][x];
                else
                    _heuristic[y][x] = 0;
            }
        }
    }

private:
    int _numAnts = 0;        // number of ants per iteration
    int _numIterations = 0;  // total iterations

    float _alpha = 0;       // pheromone importance
    float _beta = 0;        // heuristic importance
    float _rho = 0;         // evaporation rate (0–1)

    //Pheromone bounds
    float _tauMin = 0;
    float _tauMax = 0;

    int _globalBestFreq = 10;

    std::vector<Ant> _ants;
    std::vector<std::vector<float>> _pheromone;
    std::vector<std::vector<float>> _heuristic;
};
