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

    TSPMMAS(const std::vector<std::vector<int>>& adjMat, int numAnts, int numIterations, float alpha, float beta, float rho, float tauMin, float tauMax, unsigned int seed = std::random_device{}()) 
        :TSPAlgo(adjMat, seed), _numAnts(numAnts), _numIterations(numIterations), _alpha(alpha), _beta(beta), _rho(rho), _tauMin(tauMin), _tauMax(tauMax)
    {
        int nCities = adjMat.size();

        //Init pheromone matrix with _tauMax
        _pheromone.resize(nCities, std::vector<float>(nCities, _tauMax));

        //Init the heuristic matrix with 1/<dist between cities>
        _heuristic.resize(nCities, std::vector<float>(nCities));
        for (size_t y = 0; y < nCities; y++)
        {
            for (size_t x = 0; x < nCities; x++)
            {
                if (y != x && (*_adjMat)[y][x] > 0)
                    _heuristic[y][x] = 1.0f / (*_adjMat)[y][x];
                else
                    _heuristic[y][x] = 0;
            }
        }

        //Init ants
        _ants.resize(_numAnts);
    }

    void solve() override
    {

    }

private:
    int _numAnts;        // number of ants per iteration
    int _numIterations;  // total iterations

    float _alpha;       // pheromone importance
    float _beta;        // heuristic importance
    float _rho;         // evaporation rate (0–1)

    //Pheromone bounds
    float _tauMin;
    float _tauMax;

    int _globalBestFreq = 10;

    std::vector<Ant> _ants;
    std::vector<std::vector<float>> _pheromone;
    std::vector<std::vector<float>> _heuristic;
};
