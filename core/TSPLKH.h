#pragma once
#include <climits>
#include <array>
#include <cassert>
#include <deque>

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

    bool symmetrizeCandidates = true;   //Weather or not to make the candidate sets for each city symetric
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
        int to = -1;              // city this edge goes to
        long long alpha = -1;     // how promising the edge is
        long long cost = -1;      // transformed edge cost

        LKHCandidate() = default;
        LKHCandidate(int to, long long alpha, long long cost)
            : to(to), alpha(alpha), cost(cost)
        {}

        bool operator<(const LKHCandidate& other) const
        {
            if (alpha != other.alpha)
                return alpha < other.alpha;
            if (cost != other.cost)
                return cost < other.cost;
            return to < other.to;
        }
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

    //Internal representation of a city in the tour
    struct LKHTourNode
    {
        int prev = -1;                  //the previous city in the tour
        int next = -1;                  //the next city in the tour
        int rank = -1;                  //the position of the city in the tour

        LKHTourNode() = default;
        LKHTourNode(int prev, int next, int rank)
            : prev(prev), next(next), rank(rank)
        {}
    };

    //State for one sequential LK move chain.
    struct LKHMoveState
    {
        std::vector<int> t;                         //selected LK endpoints: t[0]=t1, t[1]=t2, ...
        std::vector<std::pair<int, int>> added;     //candidate edges added by the move chain
        std::vector<std::pair<int, int>> deleted;   //tour edges deleted by the move chain
        long long gain = 0;                         //current accumulated gain
        int depth = 0;                              //number of deleted tour edges in the chain
        int backtrackingCount = 0;                  //how many times this move backtracked

        void reset()
        {
            t.clear();
            added.clear();
            deleted.clear();
            gain = 0;
            depth = 0;
            backtrackingCount = 0;
        }

        void recordAddedEdge(int t1, int t2)
        {
            added.emplace_back(t1, t2);
        }

        void undoAddedEdge()
        {
            assert(!added.empty());
            added.pop_back();
        }

        void recordDeletedEdge(int t1, int t2)
        {
            deleted.emplace_back(t1, t2);
            depth++;
        }

        void undoDeletedEdge()
        {
            assert(!deleted.empty());
            assert(depth > 0);
            deleted.pop_back();
            depth--;
        }

        void pushEndpoint(int city)
        {
            t.push_back(city);
        }

        void popEndpoint()
        {
            assert(!t.empty());
            t.pop_back();
        }
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

    void addAlphaCandidate(int from, const LKHCandidate& candidate, bool bounded);

    bool isMSTEdge(int u, int v) const;

    void generateAlphaCandidates(const std::vector<std::vector<int>>& adjMat);

    void symmetrizeCandidates();

    bool isCandidateOf(int city, int candidateToCheck) const;

    //Tour manipulation funcs
    void buildInitialTourNN(const std::vector<std::vector<int>>& adjMat, int startCity = 0);

    void addCityToTour(int city, int prev, int rank);

    bool validateTourInternal(int startCity = 0) const;

    bool isEdgeInTour(int a, int b) const;

    long long calculateInternalTourCost(const std::vector<std::vector<int>>& adjMat, int startCity = 0) const;

    std::vector<int> buildOutputPath(int startCity = 0) const;  //builds the path using std::vector<int>

    void refreshTourRanks(int startCity = 0);

    void applyTwoOptMove(int t1, int t2, int t3, int t4);

    //k-opt search
    void activateAllNodes(int nCities);

    int removeFirstActiveNode();    //removes and returns the next city that should be tried in LK search

    void activateNode(int city);    //activates the city if it isnt already active    

    bool tryImproveFromNode(int t1, const std::vector<std::vector<int>>& adjMat); //This tries to find an improving LK move that starts from city t1

    bool trySequentialMove(LKHMoveState& state, const std::vector<std::vector<int>>& adjMat);

    void runLinKernighanSearch(const std::vector<std::vector<int>>& adjMat);

    bool isAddedEdge(const LKHMoveState& state, int a, int b) const; //checks weather an edge has been added in the current Search
    bool isDeletedEdge(const LKHMoveState& state, int a, int b) const; //checks weather an edge has been deleted in the current Search

    void initializeMoveState(LKHMoveState& state, int t1, int t2, long long initialGain); //It prepares the move chain after choosing the first deleted edge (t1, t2)

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

    std::vector<std::vector<LKHCandidate>> _candidates; //candidate next cities for each city by alpha nearness

    std::vector<LKHTourNode> _tourInternal;             //Internal tour representation  

    //k-opt search
    std::deque<int> _activeNodes;       //cities that still need to be tried                      
    std::vector<char> _isActive;        //flags for weather or not each city is active
};
