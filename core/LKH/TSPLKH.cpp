#include "../TSPLKH.h"

#include <chrono>
#include <iostream>

TSPLKH::TSPLKH(const LKHConfig& config, unsigned int seed)
    :TSPAlgo(seed), _config(config)
{}

void TSPLKH::solve(const std::vector<std::vector<int>>& adjMat)
{
    //num cities
    int n = adjMat.size();
    assert(n >= 3);

    //Reset values
    _bestPisSum = 0;
    _bestLowerBound = LLONG_MIN;
    _bestNorm = LLONG_MAX;

    //init containers
    _nodes.resize(n);
    _bestPis.assign(n, 0);
    _candidates.clear();
    _candidates.resize(n);
    _tourSegments.clear();
    _citySegment.assign(n, -1);
    _cityOffsetInSegment.assign(n, -1);
    _ascentCandidates.clear();
    _ascentCandidates.resize(n);
    for (int i = 0; i < n; i++)
    {
        _nodes[i].id = i;
        _nodes[i].pi = 0;
        _nodes[i].degree = 2;
        _nodes[i].lastDegree = 2;
        _nodes[i].parent = -1;
        _nodes[i].parentCost = 0;

        //ascent candidates
        _ascentCandidates[i].reserve(_config.ascentCandidates);
    }
    _piSum = 0;

    auto t0 = std::chrono::steady_clock::now();

    // Ascent and get the lower bound
    long long lowerBound = ascent(adjMat);
    auto t1 = std::chrono::steady_clock::now();

    // Build alpha-nearness candidate sets from the final ascent 1-tree.
    generateAlphaCandidates(adjMat);
    auto t2 = std::chrono::steady_clock::now();

    // Build the initial tour using NN
    buildInitialTourNN(adjMat, 0);
    auto t3 = std::chrono::steady_clock::now();

    // Run the variable-depth k-opt search.
    runLinKernighanSearch(adjMat);
    auto t4 = std::chrono::steady_clock::now();

    auto ms = [](auto a, auto b)
    {
        return std::chrono::duration_cast<std::chrono::milliseconds>(b - a).count();
    };

    std::cout << "Ascent: " << ms(t0, t1) << " ms\n";
    std::cout << "Alpha candidates: " << ms(t1, t2) << " ms\n";
    std::cout << "Initial tour: " << ms(t2, t3) << " ms\n";
    std::cout << "LK search: " << ms(t3, t4) << " ms\n";
    std::cout << "Total: " << ms(t0, t4) << " ms\n";

    //Export the tour to the output format
    _currSolution.path = buildOutputPath(0);
    _currSolution.calculateDist(adjMat);
}

long long TSPLKH::getTransformedCost(int i, int j, const std::vector<std::vector<int>>& adjMat) const
{
    return (long long)(_config.precision * (long long)adjMat[i][j] + _nodes[i].pi + _nodes[j].pi);
}

long long TSPLKH::calculateOneTreeLowerBound(long long oneTreeCost) const
{
    return oneTreeCost - (2 * _piSum);
}
