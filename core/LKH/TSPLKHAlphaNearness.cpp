#include "../TSPLKH.h"
#include <stack>

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

void TSPLKH::computeBetaValues(
    int from,
    const std::vector<std::vector<LKHTreeEdge>>& mstAdj,
    std::vector<long long>& beta) const
{
    const int n = (int)mstAdj.size();
    assert(from >= 0 && from < n);

    //Reset beta values before computing paths from this city.
    beta.assign(n, LLONG_MIN);

    //Store city, parent, and current max edge on the path from 'from'.
    struct StackState
    {
        int node;
        int parent;
        long long pathMax;
    };

    std::stack<StackState> stack;

    //The path from a city to itself has no edge, so its max is 0.
    beta[from] = 0;
    stack.push({ from, -1, 0 });

    while (!stack.empty())
    {
        const StackState state = stack.top();
        stack.pop();

        //Visit all MST neighbors, skipping the edge we came from.
        for (const auto& edge : mstAdj[state.node])
        {
            if (edge.to == state.parent)
                continue;

            //Extend the path and keep the largest edge seen so far.
            const long long nextPathMax = std::max(state.pathMax, edge.cost);
            beta[edge.to] = nextPathMax;

            stack.push({ edge.to, state.node, nextPathMax });
        }
    }

#ifndef NDEBUG
    for (int i = 0; i < n; i++)
        assert(beta[i] != LLONG_MIN);
#endif
}
