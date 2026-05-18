#include "../TSPLKH.h"
#include <algorithm>
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

    //The path from a city to itself has no edge.
    beta[from] = LLONG_MIN;
    stack.push({ from, -1, LLONG_MIN });

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
            const long long nextPathMax =
                state.pathMax == LLONG_MIN ? edge.cost : std::max(state.pathMax, edge.cost);
            beta[edge.to] = nextPathMax;

            stack.push({ edge.to, state.node, nextPathMax });
        }
    }

#ifndef NDEBUG
    for (int i = 0; i < n; i++)
        assert(i == from || beta[i] != LLONG_MIN);
#endif
}

void TSPLKH::addAlphaCandidate(
    std::vector<std::vector<LKHCandidate>>& target,
    int from,
    const LKHCandidate& candidate,
    int maxCandidates,
    bool bounded)
{
    assert(!isCandidateOf(target, from, candidate.to));

    if (bounded && maxCandidates <= 0)
        return;

    auto& candidates = target[from];

    //If the bounded list is full and this edge is not better, skip it.
    if (bounded &&
        (int)candidates.size() == maxCandidates &&
        !(candidate < candidates.back()))
        return;

    //Keep candidates sorted by alpha, then cost.
    auto pos = std::lower_bound(
        candidates.begin(),
        candidates.end(),
        candidate);

    candidates.insert(pos, candidate);

    //Keep only the best configured number of candidates.
    if (bounded)
    {
        if ((int)candidates.size() > maxCandidates)
            candidates.pop_back();
    }
}

bool TSPLKH::isMSTEdge(int u, int v) const
{
    assert(u >= 0 && u < (int)_nodes.size());
    assert(v >= 0 && v < (int)_nodes.size());

    return _nodes[u].parent == v || _nodes[v].parent == u;
}

void TSPLKH::generateAlphaCandidates(const std::vector<std::vector<int>>& adjMat)
{
    //Final LK search candidates use the smaller search candidate limit.
    generateAlphaCandidates(
        adjMat,
        _candidates,
        _config.maxCandidates,
        _config.symmetrizeCandidates);
}

void TSPLKH::generateAlphaCandidates(
    const std::vector<std::vector<int>>& adjMat,
    std::vector<std::vector<LKHCandidate>>& target,
    int maxCandidates,
    bool symmetric)
{
    const int n = (int)adjMat.size();

    assert((int)_nodes.size() == n);
    assert(n > 0);
    assert(maxCandidates > 0);

    //Reset only the target list so the same logic can build final or ascent candidates.
    target.clear();
    target.resize(n);

    //Use the current 1-tree MST to compute max-edge path values.
    const auto mstAdj = buildMSTAdjacency(adjMat);

    std::vector<long long> beta;
    const long long selectedExtraCost = getTransformedCost(_oneTreeExtraU, _oneTreeExtraV, adjMat);

    for (int from = 0; from < n; from++)
    {
        target[from].reserve(std::min(maxCandidates, n - 1));

        //beta[to] = max MST edge on the path from 'from' to 'to'.
        computeBetaValues(from, mstAdj, beta);

        for (int to = 0; to < n; to++)
        {
            if (to == from)
                continue;

            const long long cost = getTransformedCost(from, to, adjMat);
            long long alpha = 0;

            //The selected 1-tree leaf uses its special extra edge as reference.
            if (from == _oneTreeExtraU || to == _oneTreeExtraU)
            {
                alpha = isMSTEdge(from, to) ? 0 : cost - selectedExtraCost;
            }
            else
            {
                //For normal edges, alpha is the cost increase after replacing
                //the largest edge on the MST path.
                alpha = cost - beta[to];
            }

            assert(alpha >= 0);

            //Keep only the best alpha-nearest edges for this city.
            addAlphaCandidate(target, from, { to, alpha, cost }, maxCandidates, true);
        }
    }

    //Complement the candidate graph so an accepted edge is visible from both endpoints.
    if (symmetric)
        symmetrizeCandidates(target);
}

void TSPLKH::symmetrizeCandidates(std::vector<std::vector<LKHCandidate>>& target)
{
    const int n = target.size();

    // Remember original sizes so newly added reverse edges are not reprocessed.
    std::vector<int> originalSize(n);
    for (int i = 0; i < n; i++)
        originalSize[i] = target[i].size();

    // Add missing reverse edges.
    for (int from = 0; from < n; from++)
    {
        for (int idx = 0; idx < originalSize[from]; idx++)
        {
            const auto& candidate = target[from][idx];
            const int to = candidate.to;

            if (!isCandidateOf(target, to, from))
                addAlphaCandidate(target, to, { from, candidate.alpha, candidate.cost }, 0, false);
        }
    }
}

bool TSPLKH::isCandidateOf(
    const std::vector<std::vector<LKHCandidate>>& target,
    int city,
    int candidateToCheck) const
{
    const auto& candidates = target[city];

    for (const auto& candidate : candidates)
    {
        if (candidate.to == candidateToCheck)
            return true;
    }

    return false;
}
