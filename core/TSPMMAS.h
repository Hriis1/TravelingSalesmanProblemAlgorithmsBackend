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

    TSPMMAS(int numAnts, int numIterations, double alpha, double beta, double rho, unsigned int seed = std::random_device{}()) 
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

        double pBest = 0.05f;
        double pRoot = pow(pBest, 1.0f / nCities);
        double avg = nCities / 2.0f;

        _tauMin = _tauMax * (1 - pRoot) / ((avg - 1) * pRoot);

        //Init pheromone matrix with _tauMax
        _pheromone.assign(nCities, std::vector<double>(nCities, _tauMax));

        //Init the heuristic matrix with 1/<dist between cities>
        _heuristic.assign(nCities, std::vector<double>(nCities));
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
    //Ant makes a tour
    void antTour(Ant& ant, const std::vector<std::vector<int>>& adjMat)
    {
        int nCities = adjMat.size();

        while ((int)ant.path.size() < nCities)
        {
            //Choose next city
            int nextCity = chooseNextCity(ant, nCities);

            // update ant
            ant.path.push_back(nextCity);
            ant.visited[nextCity] = true;
            ant.currentCityIdx = nextCity;
        }
    }

    int chooseNextCity(Ant& ant, int nCities)
    {
        int i = ant.currentCityIdx;

        // compute weights
        double sum = 0;
        std::vector<double> weights(nCities, 0);

        for (int j = 0; j < nCities; j++) {
            if (!ant.visited[j]) {
                double tau = _pheromone[i][j];
                double eta = _heuristic[i][j];

                double w;

                if (_alpha == 1.0 && _beta == 2.0) {
                    w = tau * (eta * eta);
                }
                else {
                    w = pow(tau, _alpha) * pow(eta, _beta);
                }

                weights[j] = w;
                sum += w;
            }
        }

        // pick next city (roulette)
        double r = std::uniform_real_distribution<>(0.0, sum)(_gen);
        double cum = 0;
        int nextCity = -1;

        for (int j = 0; j < nCities; j++) {
            if (!ant.visited[j]) {
                cum += weights[j];
                if (cum >= r) {
                    nextCity = j;
                    break;
                }
            }
        }

        // safety fallback (rare floating issue)
        if (nextCity == -1) {
            for (int j = 0; j < nCities; j++) {
                if (!ant.visited[j]) {
                    nextCity = j;
                    break;
                }
            }
        }
    }

    //Normalize a path to start at the start city
    void normalizePathToStart(std::vector<int>& path, int startCity = 0)
    {
        auto it = std::find(path.begin(), path.end(), startCity);
        if (it != path.end())
            std::rotate(path.begin(), it, path.end());
    }

private:
    int _numAnts = 0;        // number of ants per iteration
    int _numIterations = 0;  // total iterations

    double _alpha = 0;       // pheromone importance
    double _beta = 0;        // heuristic importance
    double _rho = 0;         // evaporation rate (0–1)

    //Pheromone bounds
    double _tauMin = 0;
    double _tauMax = 0;

    int _globalBestFreq = 10;

    std::vector<Ant> _ants;
    std::vector<std::vector<double>> _pheromone;
    std::vector<std::vector<double>> _heuristic;
};
