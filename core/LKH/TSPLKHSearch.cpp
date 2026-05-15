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

		// this will call the variable-depth LK move search.
		// if (trySequentialMove(t1, t2, gainFromDeletedEdge, adjMat))
		//     return true;
	}

	return false;
}
