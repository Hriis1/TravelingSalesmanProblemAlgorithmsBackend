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
        int dist = 0;
        int currentCityIdx;  

        void reset(int N, int startCity) {
            dist = 0;
            path.clear();
            path.reserve(N);

            visited.assign(N, false);

            currentCityIdx = startCity;
            path.push_back(startCity);
            visited[startCity] = true;
        }
    };

public:

    TSPMMAS(int numAnts, int numIterations, double alpha, double beta, double rho, int nnoimpr = INT_MAX, unsigned int seed = std::random_device{}())
        :TSPAlgo(seed), _numAnts(numAnts), _ants(numAnts), _numIterations(numIterations), _alpha(alpha), _beta(beta), _rho(rho), _nnoimpr(nnoimpr)
    {}

    void solve(const std::vector<std::vector<int>>& adjMat) override
    {
        int nCities = adjMat.size();

        //Initial solution using nearst neighbor
        _currSolution.path = TSPUtils::nearestNeighborPath(adjMat);
        _currSolution.calculateDist(adjMat);

        //Init _tauMin and _tauMax
        double pBest = 0.05f;
        updateTauMinAndTauMax(pBest, nCities);

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

        //Iterate for _numIterations for each ant
        int startCity = 0;
        int itersWithNoImprovement = 0;
        int numItersPerGlobalBestForPheromone = 10;
        for (size_t iter = 0; iter < _numIterations; iter++)
        {
            int currIterBestDist = INT_MAX;
            int currIterBestAntIdx = -1;

            //update iters with no improvement
            ++itersWithNoImprovement;

            //Ants loop
            for (size_t a = 0; a < _numAnts; a++)
            {
                auto& currAnt = _ants[a];

                //reset ants with a random start city
                startCity = std::uniform_int_distribution<>(0, nCities - 1)(_gen);
                currAnt.reset(nCities, startCity);

                //Ant tour
                doAntTour(currAnt, adjMat);

                //2opt optimization
                twoOpt(currAnt.path, adjMat);

                //Update iteration and global best
                if (currAnt.dist < currIterBestDist)
                {
                    currIterBestDist = currAnt.dist;
                    currIterBestAntIdx = a;

                    //update global best if it beats it
                    if (currAnt.dist < _currSolution.dist)
                    {
                        //update best
                        _currSolution.path = currAnt.path;
                        _currSolution.dist = currAnt.dist;

                        //Reset iters with no improvement
                        itersWithNoImprovement = 0;
                    }
                }
            }

            //If a new global best was found
            if (itersWithNoImprovement == 0) {
                //update tauMin and tauMax
                updateTauMinAndTauMax(pBest, nCities);
            }

            //If algo got stuck with no improvements over many iterations
            if (itersWithNoImprovement >= _nnoimpr)
            {
                //reset _pheromones to _tauMax
                resetPheromones();
                itersWithNoImprovement = 0;
                continue;
            }

            //Evaporate pheromones - compares with _tauMin
            evaporatePheromone();

            //Choose iter best or global best path to deposit pheromone
            bool useGlobal = (iter % numItersPerGlobalBestForPheromone == 0);
            const std::vector<int>& depositingPath = useGlobal ? _currSolution.path : _ants[currIterBestAntIdx].path;
            int depositDist = useGlobal ? _currSolution.dist : _ants[currIterBestAntIdx].dist;

            //Deposit pheromone  - compares with _tauMax
            depositPheromone(depositingPath, depositDist);
        }
    }

private:
    //Ant makes a tour
    void doAntTour(Ant& ant, const std::vector<std::vector<int>>& adjMat)
    {
        int nCities = adjMat.size();

        while ((int)ant.path.size() < nCities)
        {
            //Choose next city
            int nextCity = chooseNextCity(ant, nCities);

            // update ant
            ant.dist += adjMat[ant.currentCityIdx][nextCity];
            ant.path.push_back(nextCity);
            ant.visited[nextCity] = true;
            ant.currentCityIdx = nextCity;
        }

        //add return to start
        ant.dist += adjMat[ant.currentCityIdx][ant.path[0]];
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

    void evaporatePheromone()
    {
        int n = _pheromone.size();
        double factor = 1.0 - _rho;

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                //evaporate
                _pheromone[i][j] *= factor;

                //clamp with the lower bound
                _pheromone[i][j] = std::max(_tauMin, _pheromone[i][j]);

                //matrix is mirrored
                _pheromone[j][i] = _pheromone[i][j];
            }
        }
    }

    void depositPheromone(const std::vector<int>& path, int dist)
    {
        //amount to deposit
        double pherAmount = 1.0 / dist;

        //deposit on path edges
        for (size_t i = 1; i < path.size(); i++)
        {
            int u = path[i - 1];
            int v = path[i];

            _pheromone[u][v] += pherAmount;

            //clamp with upper bound
            _pheromone[u][v] = std::min(_tauMax, _pheromone[u][v]);

            //mirror
            _pheromone[v][u] = _pheromone[u][v];
        }

        int u = path.back();
        int v = path[0];

        //deposit on return to start edge
        _pheromone[u][v] += pherAmount;

        //clamp with upper bound
        _pheromone[u][v] = std::min(_tauMax, _pheromone[u][v]);

        //mirror
        _pheromone[v][u] = _pheromone[u][v];
    }

    void updateTauMinAndTauMax(double pBest, int nCities)
    {
        double pRoot = pow(pBest, 1.0 / nCities);
        double avg = nCities / 2.0;

        _tauMax = 1.0 / (_rho * _currSolution.dist);
        _tauMin = _tauMax * (1 - pRoot) / ((avg - 1) * pRoot);
    }

    void resetPheromones()
    {
        int n = _pheromone.size();

        for (int i = 0; i < n; i++) {
            std::fill(_pheromone[i].begin(), _pheromone[i].end(), _tauMax);
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
    int _nnoimpr = INT_MAX;

    std::vector<Ant> _ants;
    std::vector<std::vector<double>> _pheromone;
    std::vector<std::vector<double>> _heuristic;
};
