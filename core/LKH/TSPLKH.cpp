#include "../TSPLKH.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <string>

TSPLKH::TSPLKH(const LKHConfig& config, unsigned int seed)
    :TSPAlgo(seed), _config(config)
{}

void TSPLKH::solve(const std::vector<std::vector<int>>& adjMat)
{
    //Func to display time
    auto ms = [](auto a, auto b)
    {
        return std::chrono::duration_cast<std::chrono::milliseconds>(b - a).count();
    };

    //time points
    std::chrono::steady_clock::time_point t0, t1;

    //funcs to measure time
    auto execTimeStart = [&]()
    {
        if (_config.printExecTimes)
            t0 = std::chrono::steady_clock::now();
    };
    auto execTimeEnd = [&](const std::string& label)
    {
        if (_config.printExecTimes)
        {
            t1 = std::chrono::steady_clock::now();
            std::cout << label << ": " << ms(t0, t1) << " ms\n";
        }
    };

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


    // Ascent and get the lower bound
    execTimeStart();
    long long lowerBound = ascent(adjMat);
    execTimeEnd("Ascent");


    // Build alpha-nearness candidate sets from the final ascent 1-tree.
    execTimeStart();
    generateAlphaCandidates(adjMat);
    execTimeEnd("Alpha candidates");

    //Do search trials
    const int maxTrials = std::max(1, _config.maxTrials);
    long long bestCost = LLONG_MAX;
    std::vector<int> bestPath;

    //Build the first tour once; later trials start from kicked local optima.
    execTimeStart();
    std::uniform_int_distribution<int> pickCity(0, n - 1);
    buildInitialTourWalk(0);
    execTimeEnd("Initial tour");

    for (int trial = 0; trial < maxTrials; trial++)
    {
        // Run the variable-depth k-opt search.
        execTimeStart();
        runLinKernighanSearch(adjMat);
        execTimeEnd("Search trial " + std::to_string(trial));

        //Update best cost/path
        long long cost = calculateInternalTourCost(adjMat, 0);
        if (cost < bestCost)
        {
            bestCost = cost;
            bestPath = buildOutputPath(0);
        }

        if (trial + 1 < maxTrials) //if its not last trial
        {
            if (_config.trialTourGenType == LKHTrialTourGenType::KICK) //kick generation
            {
                //Perturb the current best tour before the next trial
                execTimeStart();
                rebuildTourSegmentsFromPath(bestPath);
                applyKSwapKick(_config.kickStrength);
                execTimeEnd("Kick path - trial " + std::to_string(trial));

            }
            else //walk generation
            {
                buildInitialTourWalk(pickCity(_gen));
            }
        }
    }

    //Export the tour to the output format
    _currSolution.path = bestPath;
    _currSolution.dist = bestCost;
}

long long TSPLKH::getTransformedCost(int i, int j, const std::vector<std::vector<int>>& adjMat) const
{
    return (long long)(_config.precision * (long long)adjMat[i][j] + _nodes[i].pi + _nodes[j].pi);
}

long long TSPLKH::calculateOneTreeLowerBound(long long oneTreeCost) const
{
    return oneTreeCost - (2 * _piSum);
}
