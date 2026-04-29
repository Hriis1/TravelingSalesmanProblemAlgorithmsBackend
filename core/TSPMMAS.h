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

    TSPMMAS(int numIterations, double alpha, double beta, double rho, int nnoimpr = INT_MAX, unsigned int seed = std::random_device{}())
        :TSPAlgo(seed), _numIterations(numIterations), _alpha(alpha), _beta(beta), _rho(rho), _nnoimpr(nnoimpr)
    {}

    void solve(const std::vector<std::vector<int>>& adjMat) override
    {
        int nCities = adjMat.size();
        int nNearestNeighborsMax = 20;

        //Init
        _numAnts = nCities;
        _ants.assign(_numAnts, Ant());
        _candidates.reserve(nCities);
        _weights.reserve(nCities);

        //Compute _nearestNeighbors
        initNearestNeighborCandidates(adjMat, nNearestNeighborsMax);

        //Initial solution using nearst neighbor
        _currSolution.path = TSPUtils::nearestNeighborPath(adjMat);
        _currSolution.calculateDist(adjMat);

        //Init _tauMin and _tauMax
        double pBest = 0.05f;
        updateTauMinAndTauMax(pBest, nCities);

        //Init pheromone matrices with _tauMax
        double tauMaxPowAlpha = pow(_tauMax, _alpha);
        _pheromone.assign(nCities, std::vector<double>(nCities, _tauMax));
        _pheromonePowAlpha.assign(nCities, std::vector<double>(nCities, tauMaxPowAlpha));

        //Init the heuristi matrices with 1/<dist between cities>
        _heuristic.assign(nCities, std::vector<double>(nCities));
        _heuristicPowBeta = _heuristic;
        for (size_t y = 0; y < nCities; y++)
        {
            for (size_t x = 0; x < nCities; x++)
            {
                if (y != x && adjMat[y][x] > 0) {
                    _heuristic[y][x] = 1.0 / adjMat[y][x];
                    _heuristicPowBeta[y][x] = pow(_heuristic[y][x], _beta);
                }
                else {
                    _heuristic[y][x] = 0;
                }
            }
        }

        //Iterate for _numIterations for each ant
        int startCity = 0;
        int itersWithNoImprovement = 0;
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

            Ant& bestOfIter = _ants[currIterBestAntIdx];

            //2opt optimization for best ant or iter
            twoOpt(bestOfIter.path, adjMat);

            //Calculate the dist of the path after 2opt
            bestOfIter.dist = calculatePathDist(bestOfIter.path, adjMat);

            //update global best if it beats it
            if (bestOfIter.dist < _currSolution.dist)
            {
                //update best
                _currSolution.path = bestOfIter.path;
                _currSolution.dist = bestOfIter.dist;

                //Reset iters with no improvement
                itersWithNoImprovement = 0;
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
            bool useGlobal = (iter % _globalBestFreq == 0);
            const std::vector<int>& depositingPath = useGlobal ? _currSolution.path : _ants[currIterBestAntIdx].path;
            int depositDist = useGlobal ? _currSolution.dist : _ants[currIterBestAntIdx].dist;

            //Deposit pheromone  - compares with _tauMax
            depositPheromone(depositingPath, depositDist);
        }

        //Make city 0 the starting city of the best tour
        normalizePathToStart(_currSolution.path, 0);
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

        _candidates.clear();
        _weights.clear();

        double sum = 0;

        auto addCandidateFunc = [&](int i, int j) {
            if (!ant.visited[j]) {

                double w = _pheromonePowAlpha[i][j] * _heuristicPowBeta[i][j];

                _candidates.push_back(j);
                _weights.push_back(w);
                sum += w;
            }
        };

        //build candidates + weights from nearest neighbors
        for (size_t j = 0; j < _nearestNeighborCandidates[i].size(); j++)
        {
            int city = _nearestNeighborCandidates[i][j];
            addCandidateFunc(i, city);
        }


        //Fall back for when all nearest neighbors are visited
        if (_candidates.empty()) {
            // build candidates + weights from all cities
            for (int j = 0; j < nCities; j++) {

                addCandidateFunc(i, j);
            }
        }

        // roulette selection
        double r = std::uniform_real_distribution<>(0.0, sum)(_gen);
        double cum = 0;

        for (size_t k = 0; k < _candidates.size(); k++) {
            cum += _weights[k];
            if (cum >= r)
                return _candidates[k];
        }

        // fallback (safety)
        return _candidates.back();
    }

    //Evaporates pheromone and updates the _pheromonePowAlpha matrix
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

                //update _pheromonePowAlpha
                _pheromonePowAlpha[i][j] = pow(_pheromone[i][j], _alpha);
                _pheromonePowAlpha[j][i] = _pheromonePowAlpha[i][j];
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

            //update _pheromonePowAlpha
            _pheromonePowAlpha[u][v] = pow(_pheromone[u][v], _alpha);
            _pheromonePowAlpha[v][u] = _pheromonePowAlpha[u][v];
        }

        int u = path.back();
        int v = path[0];

        //deposit on return to start edge
        _pheromone[u][v] += pherAmount;

        //clamp with upper bound
        _pheromone[u][v] = std::min(_tauMax, _pheromone[u][v]);

        //mirror
        _pheromone[v][u] = _pheromone[u][v];

        //update _pheromonePowAlpha
        _pheromonePowAlpha[u][v] = pow(_pheromone[u][v], _alpha);
        _pheromonePowAlpha[v][u] = _pheromonePowAlpha[u][v];
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
        double tauMaxPowAlpha = pow(_tauMax, _alpha);

        for (int i = 0; i < n; i++) {
            std::fill(_pheromone[i].begin(), _pheromone[i].end(), _tauMax);
            std::fill(_pheromonePowAlpha[i].begin(), _pheromonePowAlpha[i].end(), tauMaxPowAlpha);
        }
    }

    void initNearestNeighborCandidates(const std::vector<std::vector<int>>& adjMat, int nNeighbors)
    {
        int n = adjMat.size();
        int k = std::min(nNeighbors, n - 1);

        _nearestNeighborCandidates.assign(n, std::vector<int>());

        std::vector<int> all;
        all.reserve(n - 1);

        for (int i = 0; i < n; i++)
        {
            all.clear();

            for (int j = 0; j < n; j++) {
                if (j != i)
                    all.push_back(j);
            }

            std::partial_sort(all.begin(), all.begin() + k, all.end(),
                [&](int a, int b) {
                    return adjMat[i][a] < adjMat[i][b];
                });

            _nearestNeighborCandidates[i].assign(all.begin(), all.begin() + k);
        }
    }

    //Normalize a path to start at the start city
    void normalizePathToStart(std::vector<int>& path, int startCity = 0)
    {
        auto it = std::find(path.begin(), path.end(), startCity);
        if (it != path.end())
            std::rotate(path.begin(), it, path.end());
    }

    //Unused 3opt
    void buildPos(const std::vector<int>& tour, std::vector<int>& pos)
    {
        int n = tour.size();
        pos.resize(n);

        for (int i = 0; i < n; i++)
            pos[tour[i]] = i;
    }

    bool threeOptMove(std::vector<int>& tour,
        const std::vector<std::vector<int>>& d,
        int i, int j, int k)
    {
        int n = tour.size();

        int a = tour[i];
        int b = tour[(i + 1) % n];
        int c = tour[j];
        int d1 = tour[(j + 1) % n];
        int e = tour[k];
        int f = tour[(k + 1) % n];

        int oldCost = d[a][b] + d[c][d1] + d[e][f];

        int bestDelta = 0;
        int bestCase = -1;

        int deltas[7] = {
            d[a][c] + d[b][d1] + d[e][f] - oldCost,
            d[a][b] + d[c][e] + d[d1][f] - oldCost,
            d[a][e] + d[d1][b] + d[c][f] - oldCost,
            d[a][c] + d[b][e] + d[d1][f] - oldCost,
            d[a][d1] + d[e][b] + d[c][f] - oldCost,
            d[a][e] + d[d1][c] + d[b][f] - oldCost,
            d[a][d1] + d[e][c] + d[b][f] - oldCost
        };

        for (int t = 0; t < 7; t++) {
            if (deltas[t] < bestDelta) {
                bestDelta = deltas[t];
                bestCase = t;
            }
        }

        if (bestCase == -1)
            return false;

        std::vector<int> replacement;
        replacement.reserve(k - i);

        auto appendForward = [&](int first, int last) {
            for (int idx = first; idx <= last; idx++)
                replacement.push_back(tour[idx]);
        };

        auto appendReverse = [&](int first, int last) {
            for (int idx = last; idx >= first; idx--)
                replacement.push_back(tour[idx]);
        };

        int bStart = i + 1;
        int bEnd = j;
        int cStart = j + 1;
        int cEnd = k;

        switch (bestCase)
        {
        case 0:
            appendReverse(bStart, bEnd);
            appendForward(cStart, cEnd);
            break;

        case 1:
            appendForward(bStart, bEnd);
            appendReverse(cStart, cEnd);
            break;

        case 2:
            appendReverse(cStart, cEnd);
            appendForward(bStart, bEnd);
            break;

        case 3:
            appendReverse(bStart, bEnd);
            appendReverse(cStart, cEnd);
            break;

        case 4:
            appendForward(cStart, cEnd);
            appendForward(bStart, bEnd);
            break;

        case 5:
            appendReverse(cStart, cEnd);
            appendReverse(bStart, bEnd);
            break;

        case 6:
            appendForward(cStart, cEnd);
            appendReverse(bStart, bEnd);
            break;
        }

        std::copy(replacement.begin(), replacement.end(), tour.begin() + i + 1);

        return true;
    }

    void threeOpt(std::vector<int>& tour,
        const std::vector<std::vector<int>>& adjMat)
    {
        int n = tour.size();

        std::vector<int> pos(tour.size());

        buildPos(tour, pos);

        std::vector<bool> dontLook(n, false);

        bool improved = true;

        while (improved)
        {
            improved = false;

            for (int i = 0; i < n; i++)
            {
                int a = tour[i];
                if (dontLook[a]) continue;

                bool improvedHere = false;

                for (int bCity : _nearestNeighborCandidates[a])
                {
                    int j = pos[bCity];
                    if (j <= i + 1 || j >= n - 1) continue;

                    for (int cCity : _nearestNeighborCandidates[bCity])
                    {
                        int k = pos[cCity];
                        if (k <= j + 1 || k >= n) continue;

                        if (threeOptMove(tour, adjMat, i, j, k))
                        {
                            buildPos(tour, pos);
                            std::fill(dontLook.begin(), dontLook.end(), false);

                            improved = true;
                            improvedHere = true;
                            break;
                        }
                    }

                    if (improvedHere) break;
                }

                if (!improvedHere)
                    dontLook[a] = true;
            }
        }
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
    std::vector<std::vector<double>> _pheromonePowAlpha;
    std::vector<std::vector<double>> _heuristic;
    std::vector<std::vector<double>> _heuristicPowBeta;

    //Used in chooseNextCity when ant does tours
    std::vector<std::vector<int>> _nearestNeighborCandidates;
    std::vector<int> _candidates;
    std::vector<double> _weights;
};
