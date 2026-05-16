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
		const int t2 = i == 0 ? _tourInternal[t1].prev : _tourInternal[t1].next;

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
	int t1 = state.t[0];
	int t2 = state.t[1];
	long long gain = state.gain;

	assert(isEdgeInTour(t1, t2));

	//is the edge a forward or backward move
	const bool forward = _tourInternal[t1].next == t2;

	for (const auto& candidate : _candidates[t2]) //for each of the candidate next eges of t2
	{
		const int t3 = candidate.to;

		//The added edge (t2, t3) must not already be a tour edge.
		if (isEdgeInTour(t2, t3))
			continue;

		const long long gainAfterAdd = gain - adjMat[t2][t3];
		if (gainAfterAdd <= 0)
			continue;

		const int t4 = forward ? _tourInternal[t3].prev : _tourInternal[t3].next;

		//All four endpoints must be distinct for this 2-opt close.
		if (t3 == t1 || t3 == t2 || t4 == t1 || t4 == t2)
			continue;

		const long long totalGain =
			gainAfterAdd + adjMat[t3][t4] - adjMat[t4][t1];

		if (totalGain <= 0)
			continue;

		if (forward)
			applyTwoOptMove(t1, t2, t3, t4);
		else
			applyTwoOptMove(t2, t1, t4, t3);

		activateNode(t1);
		activateNode(t2);
		activateNode(t3);
		activateNode(t4);

		return true;
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
