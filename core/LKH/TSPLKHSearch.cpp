#include "../TSPLKH.h"

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
	assert(t1 >= 0 && t1 < (int)_tourInternal.size());

	for (int i = 0; i < 2; i++)
	{
		//try edge with prev city and next city
		const int t2 = i == 0 ? getPrevInTour(t1) : getNextInTour(t1);

		if (t2 < 0)
			continue;

		//deleted edge is (t1, t2).
		const long long gainFromDeletedEdge = adjMat[t1][t2];

		//Init the LKHMoveState
		LKHMoveState moveState;
		initializeMoveState(moveState, t1, t2, gainFromDeletedEdge);

		if (trySequentialMove(moveState, adjMat))
			return true;
	}

	return false;
}

bool TSPLKH::trySequentialMove(LKHMoveState& state, const std::vector<std::vector<int>>& adjMat)
{
	assert(state.t.size() >= 2);

	int t1 = state.t[0];
	int t2 = state.t.back();
	long long initialGain = state.gain;

	if (_config.backtrackingLimit > 0 && state.backtrackingCount >= _config.backtrackingLimit)
		return false;

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
					state.recordAddedEdge(t4, t1);
					state.gain = totalGain;

					if (tryApplyKOptMove(state))
					{
						activateMoveNodes(state);
						return true;
					}

					state.gain = gainBeforeClose;
					state.undoAddedEdge();
				}
			}

			//If closing did not work, continue the sequential LK chain deeper.
			if (state.depth < _config.maxDepth && trySequentialMove(state, adjMat))
				return true;

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
	const int n = (int)_tourInternal.size();
	if (n < 3 || state.added.size() != state.deleted.size())
		return false;

	std::vector<std::array<int, 2>> adj(n);
	for (int city = 0; city < n; city++)
		adj[city] = { getPrevInTour(city), getNextInTour(city) };

	auto removeEdge = [&](int a, int b) -> bool
	{
		if (a < 0 || a >= n || b < 0 || b >= n)
			return false;

		if (adj[a][0] == b)
			adj[a][0] = -1;
		else if (adj[a][1] == b)
			adj[a][1] = -1;
		else
			return false;

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
		if (a < 0 || a >= n || b < 0 || b >= n || a == b)
			return false;

		if (adj[a][0] == b || adj[a][1] == b)
			return false;

		if (adj[b][0] == a || adj[b][1] == a)
			return false;

		int aSlot = adj[a][0] == -1 ? 0 : (adj[a][1] == -1 ? 1 : -1);
		int bSlot = adj[b][0] == -1 ? 0 : (adj[b][1] == -1 ? 1 : -1);
		if (aSlot == -1 || bSlot == -1)
			return false;

		adj[a][aSlot] = b;
		adj[b][bSlot] = a;
		return true;
	};

	//Apply the recorded edge exchange to a temporary degree-2 graph.
	for (const auto& edge : state.deleted)
	{
		if (!removeEdge(edge.first, edge.second))
			return false;
	}

	for (const auto& edge : state.added)
	{
		if (!addEdge(edge.first, edge.second))
			return false;
	}

	std::vector<LKHTourNode> newTour(n);
	std::vector<char> visited(n, 0);

	const int startCity = state.t[0];
	int prevCity = -1;
	int currCity = startCity;

	//Traverse the temporary graph; it is valid only if it forms one Hamiltonian cycle.
	for (int rank = 0; rank < n; rank++)
	{
		if (currCity < 0 || currCity >= n || visited[currCity])
			return false;

		if (adj[currCity][0] < 0 || adj[currCity][1] < 0)
			return false;

		visited[currCity] = 1;
		newTour[currCity].prev = prevCity;
		newTour[currCity].rank = rank;

		const int nextCity =
			prevCity == -1 ? adj[currCity][0] :
			adj[currCity][0] == prevCity ? adj[currCity][1] :
			adj[currCity][1] == prevCity ? adj[currCity][0] : -1;

		if (nextCity == -1)
			return false;

		newTour[currCity].next = nextCity;
		prevCity = currCity;
		currCity = nextCity;
	}

	if (currCity != startCity)
		return false;

	newTour[startCity].prev = prevCity;
	_tourInternal.swap(newTour);
	rebuildTourSegments(0);

	return validateTourInternal(0);
}

void TSPLKH::activateMoveNodes(const LKHMoveState& state)
{
	for (int city : state.t)
		activateNode(city);
}
