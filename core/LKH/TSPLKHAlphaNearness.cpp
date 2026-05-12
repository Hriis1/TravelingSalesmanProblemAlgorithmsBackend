#include "../TSPLKH.h"

std::vector<std::vector<TSPLKH::LKHTreeEdge>> TSPLKH::buildMSTAdjacency(
    const std::vector<std::vector<int>>& adjMat) const
{
    const int n = (int)adjMat.size();
    assert((int)_nodes.size() == n);

    std::vector<std::vector<LKHTreeEdge>> mstAdj(n);

    for (int i = 0; i < n; i++)
    {
        const int parent = _nodes[i].parent;
        if (parent == -1)
            continue;

        assert(parent >= 0 && parent < n);
        assert(parent != i);

        const long long cost = _nodes[i].parentCost;
        mstAdj[i].emplace_back(parent, cost);
        mstAdj[parent].emplace_back(i, cost);
    }

#ifndef NDEBUG
    int edgeCount = 0;
    int roots = 0;
    for (int i = 0; i < n; i++)
    {
        roots += _nodes[i].parent == -1 ? 1 : 0;
        edgeCount += (int)mstAdj[i].size();
    }

    assert(roots == 1);
    assert(edgeCount == 2 * (n - 1));
#endif

    return mstAdj;
}
