#include "../TSPLKH.h"

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
    for (int i = 0; i < n; i++)
    {
        _nodes[i].id = i;
        _nodes[i].pi = 0;
        _nodes[i].degree = 2;
        _nodes[i].lastDegree = 2;
        _nodes[i].parent = -1;
        _nodes[i].parentCost = 0;
    }
    _piSum = 0;

    //Ascent and get the lower bound
    long long lowerBound = ascent(adjMat);

    //Build alpha-nearness candidate sets from the final ascent 1-tree.
    generateAlphaCandidates(adjMat);

    //Build the initial tour using the internal representation using NN
    buildInitialTourNN(adjMat, 0);

    //Run the variable-depth k-opt search.
    runLinKernighanSearch(adjMat);

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
