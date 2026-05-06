#include "../TSPLKH.h"

TSPLKH::TSPLKH(const LKHConfig& config, unsigned int seed)
    :TSPAlgo(seed), _config(config)
{}

void TSPLKH::solve(const std::vector<std::vector<int>>& adjMat)
{
    //num cities
    int n = adjMat.size();

    //Reset values
    _bestPisSum = 0;
    _bestLowerBound = LLONG_MIN;
    _bestNorm = LLONG_MAX;

    //init the _nodes and reset _piSum
    _nodes.resize(n);
    _bestPis.assign(n, 0);
    for (int i = 0; i < n; i++)
    {
        _nodes[i].id = i;
        _nodes[i].pi = 0;
        _nodes[i].degree = 2;
        _nodes[i].lastDegree = 2;
        _nodes[i].parent = -1;
    }
    _piSum = 0;

    //Ascent and get the lower bound
    long long lowerBound = ascent(adjMat);
}

long long TSPLKH::getTransformedCost(int i, int j, const std::vector<std::vector<int>>& adjMat) const
{
    return (long long)(_config.precision * (long long)adjMat[i][j] + _nodes[i].pi + _nodes[j].pi);
}

long long TSPLKH::calculateOneTreeLowerBound(long long oneTreeCost) const
{
    return oneTreeCost - (2 * _piSum);
}