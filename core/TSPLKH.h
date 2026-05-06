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
    };

public:

    TSPLKH(const LKHConfig& config, unsigned int seed = std::random_device{}());

    void solve(const std::vector<std::vector<int>>& adjMat) override;


private:

    //Returns the transformed cost between 2 cities
    long long getTransformedCost(int i, int j, const std::vector<std::vector<int>>& adjMat) const;

    long long calculateOneTreeLowerBound(long long oneTreeCost) const;

    void updatePenalty(LKHNode& node, long long delta)
    {
        constexpr long long PI_MAX = INT_MAX / 10;
        constexpr long long PI_MIN = INT_MIN / 10;

        const long long oldPi = node.pi;

        node.pi += delta;

        if (node.pi > PI_MAX)
            node.pi = PI_MAX;
        else if (node.pi < PI_MIN)
            node.pi = PI_MIN;

        _piSum += node.pi - oldPi;
    }


    void saveBestPenaltyState(long long lowerBound)
    {
        assert(_bestPis.size() == _nodes.size());

        _bestLowerBound = lowerBound;
        _bestPisSum = _piSum;
        _bestNorm = _norm;

        for (size_t i = 0; i < _nodes.size(); i++)
            _bestPis[i] = _nodes[i].pi;
    }

    void restoreBestPenaltyState()
    {
        assert(_bestPis.size() == _nodes.size());

        _piSum = _bestPisSum;

        for (size_t i = 0; i < _nodes.size(); i++)
            _nodes[i].pi = _bestPis[i];
    }

    bool isOneTreeValid(int nCities) const
    {
        //A LKH-style 1-tree is an MST over all nodes plus one extra edge.
        //The MST has exactly one root with parent == -1.
        int nDegrees = 0;
        int nMstRoots = 0;
        for (int i = 0; i < nCities; i++)
        {
            nDegrees += _nodes[i].degree;
            nMstRoots += _nodes[i].parent == -1 ? 1 : 0;
        }

        //The full 1-tree has n edges, therefore total degree is 2n.
        //The MST part has n - 1 parent links, therefore one MST root.
        if (nDegrees != 2 * nCities || nMstRoots != 1)
            return false;

        return true;
    }

    std::vector<int> buildTourFromOneTree() const
    {
        const int n = (int)_nodes.size();
        assert(n >= 3);
        assert(_norm == 0);
        assert(_oneTreeExtraU != -1);
        assert(_oneTreeExtraV != -1);

        //The current 1-tree is a tour when every node has degree 2.
        //Rebuild its two-neighbor representation from the MST parent edges
        //plus the one stored non-MST extra edge.
        std::vector<std::array<int, 2>> neighbors(n, { -1, -1 });
        std::vector<int> neighborCount(n, 0);

        auto addNeighbor = [&](int u, int v)
        {
            assert(u >= 0 && u < n);
            assert(v >= 0 && v < n);
            assert(neighborCount[u] < 2);
            neighbors[u][neighborCount[u]++] = v;
        };

        auto addUndirectedEdge = [&](int u, int v)
        {
            addNeighbor(u, v);
            addNeighbor(v, u);
        };

        for (int i = 0; i < n; i++)
        {
            if (_nodes[i].parent != -1)
                addUndirectedEdge(i, _nodes[i].parent);
        }

        addUndirectedEdge(_oneTreeExtraU, _oneTreeExtraV);

        for (int i = 0; i < n; i++)
            assert(neighborCount[i] == 2);

        std::vector<int> tour;
        tour.reserve(n);

        std::vector<char> visited(n, 0);
        int prev = -1;
        int curr = 0;

        for (int step = 0; step < n; step++)
        {
            assert(!visited[curr]);
            visited[curr] = 1;
            tour.push_back(curr);

            const int a = neighbors[curr][0];
            const int b = neighbors[curr][1];
            const int next = (a != prev) ? a : b;

            prev = curr;
            curr = next;
        }

        for (int i = 0; i < n; i++)
            assert(visited[i]);
        assert(curr == tour[0]);

        return tour;
    }

    //Increases the degrees of both nodes for a selected 1-tree edge.
    void addOneTreeEdge(int u, int v)
    {
        _nodes[u].degree++;
        _nodes[v].degree++;
    }

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
    long long buildMinimumOneTree(const std::vector<std::vector<int>>& adjMat)
    {
        const int n = (int)adjMat.size();
        assert(n >= 3);
        assert((int)_nodes.size() == n);

        //Clear the previous 1-tree state and save lastdegree. Do not clear pi here
        _oneTreeExtraU = -1;
        _oneTreeExtraV = -1;
        for (int i = 0; i < n; i++)
        {
            _nodes[i].lastDegree = _nodes[i].degree;
            _nodes[i].degree = 0;
            _nodes[i].parent = -1;
        }

        long long totalCost = 0;

        //Prim's algorithm state for the MST over all nodes.
        std::vector<char> inMST(n, 0);
        std::vector<long long> bestCost(n, LLONG_MAX);
        std::vector<int> bestParent(n, -1);

        //Start Prim from node 0. It becomes the MST root and has no parent.
        bestCost[0] = 0;

        for (int step = 0; step < n; step++)
        {
            int v = -1;
            long long vCost = LLONG_MAX;

            //Pick the node outside the MST with the cheapest known
            //connection to the current MST.
            for (int i = 0; i < n; i++)
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
            //The first inserted node is the MST root and contributes no edge.
            if (bestParent[v] != -1)
            {
                const int parent = bestParent[v];
                _nodes[v].parent = parent;
                addOneTreeEdge(v, parent);
                totalCost += vCost;
            }

            //Update the cheapest known connection for every node that is
            //still outside the MST.
            for (int w = 0; w < n; w++)
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

        //Complete the 1-tree by adding one extra edge to a leaf.
        //For every MST leaf, find its cheapest non-MST extra edge.
        //LKH chooses the leaf whose such extra edge is largest.
        int selectedLeaf = -1;
        int selectedNext = -1;
        long long selectedNextCost = LLONG_MIN;

        for (int leaf = 0; leaf < n; leaf++)
        {
            if (_nodes[leaf].degree != 1)
                continue;

            //A leaf has exactly one MST neighbor. For non-root leaves it is
            //the parent. If the MST root is a leaf, its only MST neighbor is
            //its single child, so find that child and exclude it as well.
            int mstNeighbor = _nodes[leaf].parent;
            if (mstNeighbor == -1)
            {
                for (int j = 0; j < n; j++)
                {
                    if (_nodes[j].parent == leaf)
                    {
                        mstNeighbor = j;
                        break;
                    }
                }
            }
            assert(mstNeighbor != -1);

            int bestNext = -1;
            long long bestNextCost = LLONG_MAX;

            for (int j = 0; j < n; j++)
            {
                if (j == leaf || j == mstNeighbor)
                    continue;

                const long long cost = getTransformedCost(leaf, j, adjMat);
                if (cost < bestNextCost)
                {
                    bestNext = j;
                    bestNextCost = cost;
                }
            }

            assert(bestNext != -1);

            if (bestNextCost > selectedNextCost)
            {
                selectedLeaf = leaf;
                selectedNext = bestNext;
                selectedNextCost = bestNextCost;
            }
        }

        assert(selectedLeaf != -1);
        assert(selectedNext != -1);

        _oneTreeExtraU = selectedLeaf;
        _oneTreeExtraV = selectedNext;

        addOneTreeEdge(selectedLeaf, selectedNext);
        totalCost += selectedNextCost;

        //Compute norm
        _norm = 0;
        for (const auto& node : _nodes)
        {
            long long v = (long long)node.degree - 2;
            _norm += v * v;
        }

        //Check if tree is built correctly during debug
        assert(isOneTreeValid(n));

        return totalCost;
    }

    long long ascent(const std::vector<std::vector<int>>& adjMat)
    {
        //Build the initial 1-tree and get its cost
        long long oneTreeCost = buildMinimumOneTree(adjMat);
        long long lowerBound = calculateOneTreeLowerBound(oneTreeCost);

        //Set last degree = degree after 1st 1-tree
        for (auto& node : _nodes)
            node.lastDegree = node.degree;

        saveBestPenaltyState(lowerBound);

        if (_norm == 0)
        {
            _currSolution.path = buildTourFromOneTree();
            _currSolution.calculateDist(adjMat);
            return lowerBound;
        }

        //Init vars for ascent
        const int n = (int)_nodes.size();

        const int initialPeriod =
            _config.initialPeriod > 0 ? _config.initialPeriod : std::max(n / 2, 100);

        long long stepSize = _config.initialStepSize * _config.precision;
        bool initialPhase = true;

        bool stopAscent = false;

        for (int period = initialPeriod; period > 0 && stepSize > 0 && _norm != 0; period /= 2, stepSize /= 2)
        {
            for (int p = 1; stepSize > 0 && p <= period && _norm != 0; p++)
            {
                //Update penalties
                for (size_t n = 0; n < _nodes.size(); n++)
                {
                    long long v = (long long)_nodes[n].degree - 2;

                    //Only update penalty if  degree violation != 0
                    if (v != 0)
                    {
                        long long lastV = (long long)_nodes[n].lastDegree - 2;
                        long long delta = stepSize * ((7 * v) + (3 * lastV)) / 10;

                        updatePenalty(_nodes[n], delta);
                    }
                }

                //build new 1-tree with updated penalties
                oneTreeCost = buildMinimumOneTree(adjMat);
                lowerBound = calculateOneTreeLowerBound(oneTreeCost);

                //if a better lower bound was found
                if (lowerBound > _bestLowerBound ||
                    (lowerBound == _bestLowerBound && _norm < _bestNorm))
                {
                    saveBestPenaltyState(lowerBound);

                    if (initialPhase) //double stepSize on improvement during the initial phase
                        stepSize *= 2;

                    if (p == period) //double period if the improvement happens on the final iteration of the current period, capped at initialPeriod
                        period = std::min(initialPeriod, period * 2);
                }
                else  //no improvement
                {
                    // End the initial probing phase. Reset p so the reduced step size
                    // gets a full period of iterations under the normal ascent regime.
                    if (initialPhase && p > period / 2)
                    {
                        initialPhase = false;
                        p = 0;
                        stepSize = 3 * stepSize / 4;
                    }
                }

                if (_norm == 0)
                {
                    stopAscent = true;
                    break;
                }
            }

            if (stopAscent)
                break;
        }

        restoreBestPenaltyState();
        oneTreeCost = buildMinimumOneTree(adjMat);
        lowerBound = calculateOneTreeLowerBound(oneTreeCost);

        if (_norm == 0)
        {
            _currSolution.path = buildTourFromOneTree();
            _currSolution.calculateDist(adjMat);
        }

        return lowerBound;
    }

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
};
