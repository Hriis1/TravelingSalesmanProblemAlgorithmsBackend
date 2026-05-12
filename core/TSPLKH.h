#pragma once
#include <climits>
#include <array>
#include <cassert>

#include "TSPAlgo.h"
#include "TSPSolution.h"
#include "TSPUtils.h"


struct LKHConfig
{
    // ASIGNED VALUES ARE JUST DEFAULTS

    int maxTrials = 50;                 //How many attempts LKH makes
    int maxCandidates = 15;             //How many candidate edges each city considers
    int maxDepth = 5;                   //Maximum depth of the variable k-opt search
    int backtrackingLimit = 20;         //Limits how many failed alternatives are explored.
    int runs = 1;                       //How many independent full runs to do
    int kickStrength = 4;               //How strong the perturbation is when stuck

    int initialPeriod = -1;             // LKH default: max(n / 2, 100)
    int initialStepSize = 1;            // LKH default: 1

    long long precision = 100;          //precision for calculating transformed costs
};


class TSPLKH: public TSPAlgo
{
private:

    //A node is information for each city
    struct LKHNode
    {
        int id = -1;                    //Which city it is
        long long pi = 0;               //The penalty of the city              
        int degree = 0;                 //How many cities is it connected to in the 1-tree
        int lastDegree = 0;             //The degree in the prev iteration
        int parent = -1;                //Which city it connected to while building the 1-tree
        long long parentCost = 0;       //The cost to the parent
    };

    //Stores a candidate for a node
    struct LKHCandidate
    {
        int to;              // city this edge goes to
        long long alpha;     // how promising the edge is
        long long cost;      // transformed edge cost
    };

    //One edge in the MST adjacency list used for alpha-nearness path queries.
    struct LKHTreeEdge
    {
        int to = -1;
        long long cost = 0;

        LKHTreeEdge() = default;
        LKHTreeEdge(int to, long long cost)
            : to(to), cost(cost)
        {}
    };

public:

    TSPLKH(const LKHConfig& config, unsigned int seed = std::random_device{}());

    void solve(const std::vector<std::vector<int>>& adjMat) override;

private:

    //Returns the transformed cost between 2 cities
    long long getTransformedCost(int i, int j, const std::vector<std::vector<int>>& adjMat) const;

    long long calculateOneTreeLowerBound(long long oneTreeCost) const;

    void updatePenalty(LKHNode& node, long long delta);

    void saveBestPenaltyState(long long lowerBound);

    void restoreBestPenaltyState();

    bool isOneTreeValid(int nCities) const;

#ifndef NDEBUG
    bool validateMinimumOneTree(
        const std::vector<std::vector<int>>& adjMat,
        long long oneTreeCost) const;
#endif

    //Increases the degrees of both nodes for a selected 1-tree edge.
    void addOneTreeEdge(int u, int v);

    std::vector<int> buildTourFromOneTree() const;

    //Builds one minimum 1-tree using the current transformed costs.
    //
    //LKH builds a 1-tree as:
    //  1. a minimum spanning tree over all nodes
    //  2. plus one extra edge from a selected leaf
    //
    //The 1-tree structure is stored in _nodes:
    //  - _nodes[i].parent stores the MST parent
    //  - _nodes[i].degree stores the degree in the full 1-tree
    //
    //Only the transformed total cost is returned because LKH keeps the
    //tree state on the _nodes themselves.
    long long buildMinimumOneTree(const std::vector<std::vector<int>>& adjMat);

    long long ascent(const std::vector<std::vector<int>>& adjMat);

    std::vector<std::vector<LKHTreeEdge>> buildMSTAdjacency(const std::vector<std::vector<int>>& adjMat) const;

    void computeBetaValues(int from, const std::vector<std::vector<LKHTreeEdge>>& mstAdj, std::vector<long long>& beta) const;

private:
    LKHConfig _config;                  //config data for solver

    std::vector<LKHNode> _nodes;        //_nodes - each node is information about a city

    std::vector<long long> _bestPis;    //Best penalties so far for each city

    long long _piSum = 0;               //Sum of all the penalties of the nodes

    long long _bestPisSum = 0;          //The best sum of all penalties so far

    long long _bestLowerBound = LLONG_MIN; //Best lower bound so far

    long long _norm = 0;                //how far the current minimum 1-tree is from being a valid tour                
    long long _bestNorm = LLONG_MAX;    //best norm so far

    int _oneTreeExtraU = -1;            //first endpoint of the non-MST 1-tree edge
    int _oneTreeExtraV = -1;            //second endpoint of the non-MST 1-tree edge

    std::vector<std::vector<LKHCandidate>> _candidateSet; //candidate next cities for each city by alpha nearness
};
