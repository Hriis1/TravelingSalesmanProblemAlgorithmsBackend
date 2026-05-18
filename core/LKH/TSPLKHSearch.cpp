#include "../TSPLKH.h"
#include <algorithm>

void TSPLKH::activateAllNodes(int nCities)
{
	//Set each city as active
	_isActive.assign(nCities, 1);

	//Add them to the queue
	_activeNodes.clear();
	for (size_t i = 0; i < nCities; i++)
		_activeNodes.push_back(i);
}

int TSPLKH::removeFirstActiveNode()
{
	if (_activeNodes.empty())
		return -1;

	int city = _activeNodes.front();
	_activeNodes.pop_front();
	_isActive[city] = 0;

	return city;
}

void TSPLKH::activateNode(int city)
{
	if (_isActive[city])
		return;

	_isActive[city] = 1;
	_activeNodes.push_back(city);
}

bool TSPLKH::tryImproveFromNode(int t1, const std::vector<std::vector<int>>& adjMat)
{
	assert(t1 >= 0 && t1 < (int)_citySegment.size());

	for (int i = 0; i < 2; i++)
	{
		//try edge with prev city and next city
		const int t2 = i == 0 ? getPrevInTour(t1) : getNextInTour(t1);

		if (t2 < 0)
			continue;

		//A BestKOptMove attempt may apply temporary non-improving moves.
		//Keep the original tour so the whole attempt can be restored if it fails.
		const std::vector<int> originalPath = buildOutputPath(0);

		int currentT2 = t2;
		long long currentGain = adjMat[t1][t2];
		std::vector<std::pair<int, int>> excludedEdges;
		std::vector<LKHMoveState> appliedPromisingMoves;

		for (int swaps = 0; swaps < (int)_citySegment.size(); swaps++)
		{
			LKHMoveState moveState;
			initializeMoveState(moveState, t1, currentT2, currentGain);

			LKHBestMoveResult result;
			findBestKOptMove(moveState, result, adjMat, excludedEdges);

			if (result.foundImprovement)
			{
				//A positive feasible close finishes the whole LK chain.
				if (!tryApplyKOptMove(result.improvingMove))
					break;

				for (const auto& move : appliedPromisingMoves)
					activateMoveNodes(move);
				activateMoveNodes(result.improvingMove);
				return true;
			}

			if (!result.foundPromisingMove)
				break;

			//Apply the best feasible non-improving K-opt move and continue
			//by deleting its closing edge in the next BestKOptMove call.
			if (!tryApplyKOptMove(result.promisingMove))
				break;

			appliedPromisingMoves.push_back(result.promisingMove);

			//LKH excludes deleted edges from the temporary chain to avoid
			//cycling back through the same exchange.
			for (const auto& edge : result.promisingMove.deleted)
				excludedEdges.push_back(edge);

			currentT2 = result.nextEndpoint;
			currentGain = result.promisingGain;
		}

		//No final improvement was found from this first edge.
		rebuildTourSegmentsFromPath(originalPath);
	}

	return false;
}

bool TSPLKH::findBestKOptMove(
	LKHMoveState& state,
	LKHBestMoveResult& result,
	const std::vector<std::vector<int>>& adjMat,
	const std::vector<std::pair<int, int>>& excludedEdges)
{
	assert(state.t.size() >= 2);

	int t1 = state.t[0];
	int t2 = state.t.back();
	long long initialGain = state.gain;

	if (_config.backtrackingLimit > 0 && state.backtrackingCount >= _config.backtrackingLimit)
		return false;

	auto isExcluded = [&](int a, int b)
	{
		for (const auto& edge : excludedEdges)
		{
			if (edge.first == a && edge.second == b)
				return true;
			if (edge.first == b && edge.second == a)
				return true;
		}

		return false;
	};

	for (const auto& candidate : _candidates[t2]) //for each of the candidate next eges of t2
	{
		const int t3 = candidate.to;

		//The added edge (t2, t3) must not already be a tour edge.
		if (isEdgeInTour(t2, t3))
			continue;

		//A city can only appear once as a move endpoint before the final close to t1.
		if (isEndpointUsed(state, t3))
			continue;

		const long long gainAfterAdd = initialGain - adjMat[t2][t3];
		if (gainAfterAdd <= 0)
			continue;

		//Record the added edge so deeper LK steps know this branch already uses it.
		state.recordAddedEdge(t2, t3);
		state.pushEndpoint(t3);
		state.gain = gainAfterAdd;

		const int tourNeighbors[2] = { getPrevInTour(t3), getNextInTour(t3) };

		for (int i = 0; i < 2; i++)
		{
			const int t4 = tourNeighbors[i];

			//The next deleted edge must lead to a new endpoint.
			if (isEndpointUsed(state, t4))
				continue;

			//Do not reuse edges excluded by earlier temporary K-opt moves.
			if (isExcluded(t3, t4))
				continue;

			//Record the deleted tour edge so the branch now ends at t4.
			state.recordDeletedEdge(t3, t4);
			state.pushEndpoint(t4);
			state.gain += adjMat[t3][t4];

			const long long gainBeforeClose = state.gain;

			//Try closing the chain by adding the edge from the current endpoint to t1.
			if (!isEdgeInTour(t4, t1) && !isAddedEdge(state, t4, t1))
			{
				const long long totalGain = state.gain - adjMat[t4][t1];
				if (totalGain > 0)
				{
					//Temporarily record the closing edge so this state is a complete k-opt move.
					state.recordAddedEdge(t4, t1);
					state.gain = totalGain;

					//LKH applies the first feasible positive close immediately.
					std::vector<int> unusedPath;
					if (buildPathAfterKOptMove(state, unusedPath))
					{
						result.foundImprovement = true;
						result.improvingMove = state;
						result.improvementGain = totalGain;

						state.gain = gainBeforeClose;
						state.undoAddedEdge();
						state.popEndpoint();
						state.undoDeletedEdge();
						state.gain = gainAfterAdd;
						state.popEndpoint();
						state.undoAddedEdge();
						return true;
					}

					//Undo the closing edge before continuing this branch deeper.
					state.gain = gainBeforeClose;
					state.undoAddedEdge();
				}
			}

			//Continue the sequential LK chain deeper to look for an even better close.
			if (state.depth < _config.maxDepth)
			{
				if (findBestKOptMove(state, result, adjMat, excludedEdges))
				{
					state.gain = gainAfterAdd;
					state.popEndpoint();
					state.undoDeletedEdge();
					state.gain = initialGain;
					state.popEndpoint();
					state.undoAddedEdge();
					return true;
				}
			}
			else if (t4 != t1 && !isAddedEdge(state, t4, t1) && !isEdgeInTour(t4, t1))
			{
				const long long totalGain = gainBeforeClose - adjMat[t4][t1];

				//At maximum depth, LKH remembers the best feasible non-improving
				//K-opt move and uses it as a continuation if no improvement exists.
				if (totalGain <= 0 && gainBeforeClose > result.promisingGain)
				{
					state.recordAddedEdge(t4, t1);
					state.gain = totalGain;

					std::vector<int> unusedPath;
					if (buildPathAfterKOptMove(state, unusedPath))
					{
						result.foundPromisingMove = true;
						result.promisingMove = state;
						result.promisingGain = gainBeforeClose;
						result.nextEndpoint = t4;
					}

					state.gain = gainBeforeClose;
					state.undoAddedEdge();
				}
			}

			//This branch failed, so undo the last deleted edge before trying the next one.
			state.backtrackingCount++;
			state.gain = gainAfterAdd;
			state.popEndpoint();
			state.undoDeletedEdge();

			if (_config.backtrackingLimit > 0 && state.backtrackingCount >= _config.backtrackingLimit)
				break;
		}

		//Undo the added edge before trying the next candidate.
		state.gain = initialGain;
		state.popEndpoint();
		state.undoAddedEdge();

		if (_config.backtrackingLimit > 0 && state.backtrackingCount >= _config.backtrackingLimit)
			return false;
	}

	return false;
}

void TSPLKH::runLinKernighanSearch(const std::vector<std::vector<int>>& adjMat)
{
	int n = adjMat.size();

	activateAllNodes(n); //activate all cities

	//While there are active cities try an improving move on them
	while (true)
	{
		int t1 = removeFirstActiveNode();
		if (t1 == -1)
			break;

		tryImproveFromNode(t1, adjMat);
	}
}

bool TSPLKH::isAddedEdge(const LKHMoveState& state, int a, int b) const
{
	for (size_t i = 0; i < state.added.size(); i++)
	{
		//check both (a, b) edge and (b, a) edge
		if (state.added[i].first == a && state.added[i].second == b)
			return true;

		if (state.added[i].second == a && state.added[i].first == b)
			return true;
	}

	return false;
}

bool TSPLKH::isDeletedEdge(const LKHMoveState& state, int a, int b) const
{
	for (size_t i = 0; i < state.deleted.size(); i++)
	{
		//check both (a, b) edge and (b, a) edge
		if (state.deleted[i].first == a && state.deleted[i].second == b)
			return true;

		if (state.deleted[i].second == a && state.deleted[i].first == b)
			return true;
	}

	return false;
}

void TSPLKH::initializeMoveState(LKHMoveState& state, int t1, int t2, long long initialGain)
{
	state.reset();
	state.pushEndpoint(t1);
	state.pushEndpoint(t2);
	state.recordDeletedEdge(t1, t2);
	state.gain += initialGain;
}

bool TSPLKH::isEndpointUsed(const LKHMoveState& state, int city) const
{
	for (int endpoint : state.t)
	{
		if (endpoint == city)
			return true;
	}

	return false;
}

bool TSPLKH::tryApplyKOptMove(const LKHMoveState& state)
{
	std::vector<int> newPath;
	if (!buildPathAfterKOptMove(state, newPath))
		return false;

	//Make the segment tour the primary accepted representation.
	rebuildTourSegmentsFromPath(newPath);

	return validateTourStructure(0);
}

bool TSPLKH::buildPathAfterKOptMove(const LKHMoveState& state, std::vector<int>& newPath) const
{
	const int n = (int)_citySegment.size();
	if (n < 3 || state.added.size() != state.deleted.size())
		return false;

	//Build a temporary degree-2 graph from the current segment-based tour.
	std::vector<std::array<int, 2>> adj(n);
	for (int city = 0; city < n; city++)
		adj[city] = { getPrevInTour(city), getNextInTour(city) };

	auto removeEdge = [&](int a, int b) -> bool
	{
		//Both endpoints must be valid cities.
		if (a < 0 || a >= n || b < 0 || b >= n)
			return false;

		//Remove b from a's two graph neighbors.
		if (adj[a][0] == b)
			adj[a][0] = -1;
		else if (adj[a][1] == b)
			adj[a][1] = -1;
		else
			return false;

		//Remove a from b's two graph neighbors.
		if (adj[b][0] == a)
			adj[b][0] = -1;
		else if (adj[b][1] == a)
			adj[b][1] = -1;
		else
			return false;

		return true;
	};

	auto addEdge = [&](int a, int b) -> bool
	{
		//Reject invalid cities and self-edges.
		if (a < 0 || a >= n || b < 0 || b >= n || a == b)
			return false;

		//The edge must not already exist in the temporary graph.
		if (adj[a][0] == b || adj[a][1] == b)
			return false;

		if (adj[b][0] == a || adj[b][1] == a)
			return false;

		//Each endpoint must have one free degree slot after deletions.
		int aSlot = adj[a][0] == -1 ? 0 : (adj[a][1] == -1 ? 1 : -1);
		int bSlot = adj[b][0] == -1 ? 0 : (adj[b][1] == -1 ? 1 : -1);
		if (aSlot == -1 || bSlot == -1)
			return false;

		adj[a][aSlot] = b;
		adj[b][bSlot] = a;
		return true;
	};

	//First delete the tour edges selected by the LK move.
	for (const auto& edge : state.deleted)
	{
		if (!removeEdge(edge.first, edge.second))
			return false;
	}

	//Then add the candidate/closing edges selected by the LK move.
	for (const auto& edge : state.added)
	{
		if (!addEdge(edge.first, edge.second))
			return false;
	}

	newPath.clear();
	newPath.reserve(n);

	std::vector<char> visited(n, 0);

	const int startCity = state.t[0];
	int prevCity = -1;
	int currCity = startCity;

	//Traverse the temporary graph and collect the accepted tour path.
	for (int step = 0; step < n; step++)
	{
		if (currCity < 0 || currCity >= n || visited[currCity])
			return false;

		if (adj[currCity][0] < 0 || adj[currCity][1] < 0)
			return false;

		visited[currCity] = 1;
		newPath.push_back(currCity);

		//Continue through the neighbor that is not the city we came from.
		const int nextCity =
			prevCity == -1 ? adj[currCity][0] :
			adj[currCity][0] == prevCity ? adj[currCity][1] :
			adj[currCity][1] == prevCity ? adj[currCity][0] : -1;

		if (nextCity == -1)
			return false;

		prevCity = currCity;
		currCity = nextCity;
	}

	if (currCity != startCity)
		return false;

	return true;
}

void TSPLKH::activateMoveNodes(const LKHMoveState& state)
{
	for (int city : state.t)
		activateNode(city);
}
