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
