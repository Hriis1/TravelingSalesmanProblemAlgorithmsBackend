#pragma once
#include <climits>
#include <array>
#include <cassert>
#include <queue>

#include "TSPAlgo.h"
#include "TSPSolution.h"
#include "TSPUtils.h"


struct LKHConfig
{
    // ASIGNED VALUES ARE JUST DEFAULTS

    int maxTrials = 50;                 //How many attempts LKH makes
    int maxCandidates = 15;             //How many candidate edges each city considers
    int maxDepth = 3;                   //Maximum depth of the variable k-opt search
    int backtrackingLimit = 20;         //Limits how many failed alternatives are explored.
    int runs = 1;                       //How many independent full runs to do
    int kickStrength = 4;               //How strong the perturbation is when stuck


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
        int parent = -1;                //Which city it connected to while building the 1-tree
    };

public:

    TSPLKH(const LKHConfig& config, unsigned int seed = std::random_device{}())
        :TSPAlgo(seed), _config(config)
    {}

    void solve(const std::vector<std::vector<int>>& adjMat) override
    {
        //num cities
        int n = adjMat.size();

        //init the _nodes and reset _piSum
        _nodes.resize(n);
        for (int i = 0; i < n; i++)
        {
            _nodes[i].id = i;
            _nodes[i].pi = 0;
            _nodes[i].degree = 0;
            _nodes[i].parent = -1;
        }
        _piSum = 0;

        //Build the initial 1-tree and get its cost
        long long oneTreeCost = buildMinimumOneTree(adjMat);
    }

private:

    //Returns the transformed cost between 2 cities
    long long getTransformedCost(int i, int j, const std::vector<std::vector<int>>& adjMat) const
    {
        return (long long)(_config.precision * (long long)adjMat[i][j] + _nodes[i].pi + _nodes[j].pi);
    }

    long long calculateOneTreeLowerBound(long long oneTreeCost) const
    {
        return oneTreeCost - (2 * _piSum);
    }

    void updatePenalty(LKHNode& node, long long delta)
    {
        node.pi += delta;
        _piSum += delta;
    }

    //Builds one minimum 1-tree using the current transformed costs.
    //
    //A 1-tree is:
    //  1. a minimum spanning tree over all non-root _nodes
    //  2. plus the two cheapest edges from the root into that tree
    //
    //The 1-tree structure is stored in _nodes:
    //  - _nodes[i].parent stores the MST parent for non-root _nodes
    //  - _nodes[i].degree stores the degree in the full 1-tree
    //
    //Only the transformed total cost is returned because LKH keeps the
    //tree state on the _nodes themselves.
    long long buildMinimumOneTree(const std::vector<std::vector<int>>& adjMat)
    {
        const int n = (int)adjMat.size();
        assert(n >= 3);
        assert((int)_nodes.size() == n);

        //LKH builds the 1-tree by excluding one special root from the MST.
        const int root = 0;

        //Clear the previous 1-tree state. Do not clear pi here
        for (int i = 0; i < n; i++)
        {
            _nodes[i].degree = 0;
            _nodes[i].parent = -1;
        }

        long long totalCost = 0;

        //Prim's algorithm state for the MST over _nodes 1..n-1.
        //The root is excluded from this MST and gets attached afterward.
        std::vector<char> inMST(n, 0);
        std::vector<long long> bestCost(n, LLONG_MAX);
        std::vector<int> bestParent(n, -1);

        //Start Prim from node 1. It enters the MST with no incoming edge.
        bestCost[1] = 0;

        for (int step = 1; step < n; step++)
        {
            int v = -1;
            long long vCost = LLONG_MAX;

            //Pick the non-root node outside the MST with the cheapest
            //known connection to the current MST.
            for (int i = 1; i < n; i++)
            {
                if (!inMST[i] && bestCost[i] < vCost)
                {
                    v = i;
                    vCost = bestCost[i];
                }
            }

            assert(v != -1);
            inMST[v] = 1;

            //If v has a parent, add that edge to the MST part of the 1-tree.
            //The first inserted node has no parent and contributes no edge.
            if (bestParent[v] != -1)
            {
                const int parent = bestParent[v];
                _nodes[v].parent = parent;
                _nodes[v].degree++;
                _nodes[parent].degree++;
                totalCost += vCost;
            }

            //Update the cheapest known connection for every non-root node
            //that is still outside the MST.
            for (int w = 1; w < n; w++)
            {
                if (inMST[w] || w == v)
                    continue;

                const long long cost = getTransformedCost(v, w, adjMat);
                if (cost < bestCost[w])
                {
                    bestCost[w] = cost;
                    bestParent[w] = v;
                }
            }
        }

        //Complete the 1-tree by adding the two cheapest transformed-cost
        //edges from the excluded root to two distinct non-root _nodes.
        int first = -1;
        int second = -1;
        long long firstCost = LLONG_MAX;
        long long secondCost = LLONG_MAX;

        for (int i = 1; i < n; i++)
        {
            const long long cost = getTransformedCost(root, i, adjMat);

            if (cost < firstCost)
            {
                second = first;
                secondCost = firstCost;
                first = i;
                firstCost = cost;
            }
            else if (cost < secondCost)
            {
                second = i;
                secondCost = cost;
            }
        }

        assert(first != -1);
        assert(second != -1);
        assert(first != second);

        _nodes[root].degree += 2;
        _nodes[first].degree++;
        _nodes[second].degree++;
        totalCost += firstCost + secondCost;

        return totalCost;
    }

private:
    LKHConfig _config; //config data for solver

    std::vector<LKHNode> _nodes; //_nodes - each node is information about a city

    long long _piSum = 0; //Sum of all the penalties of the nodes



};
