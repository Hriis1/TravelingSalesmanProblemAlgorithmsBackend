#include "../TSPLKH.h"
#include <queue>

void TSPLKH::updatePenalty(LKHNode& node, long long delta)
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

void TSPLKH::saveBestPenaltyState(long long lowerBound)
{
    assert(_bestPis.size() == _nodes.size());

    _bestLowerBound = lowerBound;
    _bestPisSum = _piSum;
    _bestNorm = _norm;

    for (size_t i = 0; i < _nodes.size(); i++)
        _bestPis[i] = _nodes[i].pi;
}

void TSPLKH::restoreBestPenaltyState()
{
    assert(_bestPis.size() == _nodes.size());

    _piSum = _bestPisSum;

    for (size_t i = 0; i < _nodes.size(); i++)
        _nodes[i].pi = _bestPis[i];
}

bool TSPLKH::isOneTreeValid(int nCities) const
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

#ifndef NDEBUG
bool TSPLKH::validateMinimumOneTree(
    const std::vector<std::vector<int>>& adjMat,
    long long oneTreeCost) const
{
    const int n = (int)adjMat.size();
    assert(n >= 3);
    assert((int)_nodes.size() == n);
    assert(_oneTreeExtraU >= 0 && _oneTreeExtraU < n);
    assert(_oneTreeExtraV >= 0 && _oneTreeExtraV < n);
    assert(_oneTreeExtraU != _oneTreeExtraV);
    assert(isOneTreeValid(n));

    std::vector<std::vector<char>> inMST(n, std::vector<char>(n, 0));
    std::vector<int> mstDegree(n, 0);
    long long recomputedCost = 0;
    int roots = 0;
    int mstEdges = 0;

    for (int i = 0; i < n; i++)
    {
        const int parent = _nodes[i].parent;
        if (parent == -1)
        {
            roots++;
            continue;
        }

        assert(parent >= 0 && parent < n);
        assert(parent != i);
        assert(!inMST[i][parent]);

        inMST[i][parent] = 1;
        inMST[parent][i] = 1;
        mstDegree[i]++;
        mstDegree[parent]++;
        recomputedCost += getTransformedCost(i, parent, adjMat);
        mstEdges++;
    }

    assert(roots == 1);
    assert(mstEdges == n - 1);
    assert(mstDegree[_oneTreeExtraU] == 1);
    assert(!inMST[_oneTreeExtraU][_oneTreeExtraV]);

    recomputedCost += getTransformedCost(_oneTreeExtraU, _oneTreeExtraV, adjMat);
    assert(recomputedCost == oneTreeCost);

    long long selectedLeafBestExtra = LLONG_MAX;
    for (int j = 0; j < n; j++)
    {
        if (j == _oneTreeExtraU || inMST[_oneTreeExtraU][j])
            continue;

        selectedLeafBestExtra = std::min(
            selectedLeafBestExtra,
            getTransformedCost(_oneTreeExtraU, j, adjMat));
    }
    assert(selectedLeafBestExtra == getTransformedCost(_oneTreeExtraU, _oneTreeExtraV, adjMat));

    for (int leaf = 0; leaf < n; leaf++)
    {
        if (mstDegree[leaf] != 1)
            continue;

        long long leafBestExtra = LLONG_MAX;
        for (int j = 0; j < n; j++)
        {
            if (j == leaf || inMST[leaf][j])
                continue;

            leafBestExtra = std::min(
                leafBestExtra,
                getTransformedCost(leaf, j, adjMat));
        }

        assert(leafBestExtra != LLONG_MAX);
        assert(leafBestExtra <= selectedLeafBestExtra);
    }

    long long recomputedNorm = 0;
    for (int i = 0; i < n; i++)
    {
        const long long v = (long long)_nodes[i].degree - 2;
        recomputedNorm += v * v;
        assert(_nodes[i].degree == mstDegree[i] +
            (i == _oneTreeExtraU ? 1 : 0) +
            (i == _oneTreeExtraV ? 1 : 0));
    }
    assert(recomputedNorm == _norm);

    return true;
}
#endif

void TSPLKH::addOneTreeEdge(int u, int v)
{
    _nodes[u].degree++;
    _nodes[v].degree++;
}

std::vector<int> TSPLKH::buildTourFromOneTree() const
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

long long TSPLKH::buildMinimumOneTree(const std::vector<std::vector<int>>& adjMat, bool sparse)
{
    const int n = (int)adjMat.size();
    assert(n >= 3);
    assert((int)_nodes.size() == n);

    const bool useSparse = sparse && !_ascentCandidates.empty();

    //Save previous degrees once. If sparse construction falls back to dense,
    //lastDegree must still refer to the tree from the previous ascent step.
    std::vector<int> previousDegree(n);
    for (int i = 0; i < n; i++)
        previousDegree[i] = _nodes[i].degree;

    auto resetOneTreeState = [&]()
    {
        _oneTreeExtraU = -1;
        _oneTreeExtraV = -1;
        for (int i = 0; i < n; i++)
        {
            _nodes[i].lastDegree = previousDegree[i];
            _nodes[i].degree = 0;
            _nodes[i].parent = -1;
            _nodes[i].parentCost = 0;
        }
    };

    resetOneTreeState();

    //Prim's algorithm state for the MST over all nodes.
    std::vector<char> inMST(n, 0);
    std::vector<long long> bestCost(n, LLONG_MAX);
    std::vector<int> bestParent(n, -1);

    long long totalCost = 0;
    int nodesInMST = 0;
    bool sparseTreeBuilt = false;

    if (useSparse)
    {
        using QueueItem = std::pair<long long, int>;
        std::priority_queue<QueueItem, std::vector<QueueItem>, std::greater<QueueItem>> queue;

        //Sparse Prim starts from city 0 and only relaxes ascent-candidate edges.
        bestCost[0] = 0;
        queue.emplace(0, 0);

        while (!queue.empty() && nodesInMST < n)
        {
            const QueueItem top = queue.top();
            const long long vCost = top.first;
            const int v = top.second;
            queue.pop();

            //Skip stale heap entries that were improved after insertion.
            if (inMST[v] || vCost != bestCost[v])
                continue;

            inMST[v] = 1;
            nodesInMST++;

            //If v has a parent, add that candidate edge to the sparse MST.
            if (bestParent[v] != -1)
            {
                const int parent = bestParent[v];
                _nodes[v].parent = parent;
                _nodes[v].parentCost = vCost;
                addOneTreeEdge(v, parent);
                totalCost += vCost;
            }

            //Only inspect ascent candidates. This is the LKH-style speedup.
            for (const auto& candidate : _ascentCandidates[v])
            {
                const int w = candidate.to;
                if (w < 0 || w >= n || inMST[w])
                    continue;

                const long long cost = getTransformedCost(v, w, adjMat);
                if (cost < bestCost[w])
                {
                    bestCost[w] = cost;
                    bestParent[w] = v;
                    queue.emplace(cost, w);
                }
            }
        }

        //If the sparse graph is disconnected, rebuild this 1-tree densely.
        //That keeps the state valid; the LKH-style quality safeguard comes next.
        if (nodesInMST == n)
        {
            sparseTreeBuilt = true;
        }
        else
        {
            resetOneTreeState();
            inMST.assign(n, 0);
            bestCost.assign(n, LLONG_MAX);
            bestParent.assign(n, -1);
            totalCost = 0;
            nodesInMST = 0;
        }
    }

    if (!useSparse || nodesInMST != n)
    {
        //Start dense Prim from city 0. It becomes the MST root and has no parent.
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
            nodesInMST++;

            //If v has a parent, add that edge to the MST part of the 1-tree.
            //The first inserted node is the MST root and contributes no edge.
            if (bestParent[v] != -1)
            {
                const int parent = bestParent[v];
                _nodes[v].parent = parent;
                _nodes[v].parentCost = vCost;
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

        if (sparseTreeBuilt)
        {
            //In sparse ascent mode, the extra edge is also selected from
            //the ascent candidate graph, matching the sparse 1-tree idea.
            for (const auto& candidate : _ascentCandidates[leaf])
            {
                const int j = candidate.to;
                if (j == leaf || j == mstNeighbor)
                    continue;

                const long long cost = getTransformedCost(leaf, j, adjMat);
                if (cost < bestNextCost)
                {
                    bestNext = j;
                    bestNextCost = cost;
                }
            }
        }
        else
        {
            //Dense mode scans the complete graph for the cheapest valid extra edge.
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
        }

        //If a sparse leaf lacks a usable non-MST candidate, scan all edges
        //for that leaf so the 1-tree remains valid.
        if (bestNext == -1 && sparseTreeBuilt)
        {
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

#ifndef NDEBUG
    if (!sparseTreeBuilt)
        assert(validateMinimumOneTree(adjMat, totalCost));
#endif

    return totalCost;
}

long long TSPLKH::ascent(const std::vector<std::vector<int>>& adjMat)
{
    const int n = (int)_nodes.size();
    int ascentCandidateCount = std::min(_config.ascentCandidates, n);

    const int initialPeriod =
        _config.initialPeriod > 0 ? _config.initialPeriod : std::max(n / 2, 100);

    while (true)
    {
        //A sparse-graph restart follows LKH: reset Pi and redo ascent
        //with a larger ascent candidate graph.
        _piSum = 0;
        _bestPisSum = 0;
        _bestLowerBound = LLONG_MIN;
        _bestNorm = LLONG_MAX;
        for (auto& node : _nodes)
            node.pi = 0;

        //Build the first 1-tree in the full graph.
        long long oneTreeCost = buildMinimumOneTree(adjMat);
        long long lowerBound = calculateOneTreeLowerBound(oneTreeCost);
        long long safeguardLowerBound = lowerBound;

        if (_norm == 0)
        {
            _currSolution.path = buildTourFromOneTree();
            _currSolution.calculateDist(adjMat);
            return lowerBound;
        }

        saveBestPenaltyState(lowerBound);

        //Generate the sparse graph used by repeated ascent 1-tree builds.
        if (ascentCandidateCount > 0)
            generateAlphaCandidates(adjMat, _ascentCandidates, ascentCandidateCount, true);
        else
            _ascentCandidates.clear();

        //Use this first 1-tree as the previous degree vector for the first update.
        for (auto& node : _nodes)
            node.lastDegree = node.degree;

        long long stepSize = _config.initialStepSize * _config.precision;
        bool initialPhase = true;
        bool stopAscent = false;
        bool restartAscent = false;

        for (int period = initialPeriod; period > 0 && stepSize > 0 && _norm != 0; period /= 2, stepSize /= 2)
        {
            for (int p = 1; stepSize > 0 && p <= period && _norm != 0; p++)
            {
                //Update penalties in the Held-Karp subgradient direction.
                for (size_t i = 0; i < _nodes.size(); i++)
                {
                    long long v = (long long)_nodes[i].degree - 2;

                    //Only cities with degree violation change their Pi value.
                    if (v != 0)
                    {
                        long long lastV = (long long)_nodes[i].lastDegree - 2;
                        long long delta = stepSize * ((7 * v) + (3 * lastV)) / 10;

                        updatePenalty(_nodes[i], delta);
                    }
                }

                //Most ascent iterations use the sparse ascent-candidate graph.
                oneTreeCost = buildMinimumOneTree(adjMat, ascentCandidateCount > 0);
                lowerBound = calculateOneTreeLowerBound(oneTreeCost);

                bool improved =
                    lowerBound > _bestLowerBound ||
                    (lowerBound == _bestLowerBound && _norm < _bestNorm);

                //LKH safeguard: if a sparse bound looks unrealistically high,
                //verify it with a full 1-tree before trusting it.
                if (improved &&
                    ascentCandidateCount > 0 &&
                    ascentCandidateCount < n &&
                    lowerBound > 2 * safeguardLowerBound)
                {
                    oneTreeCost = buildMinimumOneTree(adjMat, false);
                    lowerBound = calculateOneTreeLowerBound(oneTreeCost);

                    //If the full graph gives a worse lower bound than the old
                    //safeguard baseline, the sparse graph was too small.
                    if (lowerBound < safeguardLowerBound)
                    {
                        ascentCandidateCount = std::min(n, 2 * ascentCandidateCount);
                        restartAscent = true;
                        break;
                    }

                    safeguardLowerBound = lowerBound;
                    improved =
                        lowerBound > _bestLowerBound ||
                        (lowerBound == _bestLowerBound && _norm < _bestNorm);
                }

                if (improved)
                {
                    saveBestPenaltyState(lowerBound);

                    //During the initial probing phase, LKH grows the step faster
                    //when the bound keeps improving.
                    if (initialPhase)
                        stepSize *= 2;

                    //If improvement happens at the end of the period, try a longer period.
                    if (p == period)
                        period = std::min(initialPeriod, period * 2);
                }
                else
                {
                    //End the initial probing phase and continue with a reduced step.
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

            if (restartAscent || stopAscent)
                break;
        }

        if (restartAscent)
            continue;

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
}
