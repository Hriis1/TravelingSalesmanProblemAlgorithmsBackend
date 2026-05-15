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

		if (trySequentialMove(t1, t2, gainFromDeletedEdge, adjMat))
			return true;
	}

	return false;
}

bool TSPLKH::trySequentialMove(int t1, int t2, long long gain, const std::vector<std::vector<int>>& adjMat)
{
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
