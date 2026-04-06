#pragma once
#include <random>
#include <vector>

#include "TSPSolution.h"

class TSPAlgo
{
public:
	TSPAlgo(const std::vector<std::vector<int>>& adjMat, unsigned int seed)
		: _adjMat(&adjMat), _currSolution(adjMat.size()), _gen(seed)
	{}

	void reseed(unsigned int seed = std::random_device{}())
	{
		_gen.seed(seed);
	}

	virtual void solve() = 0;

	//getters
	int getCurrSolutionDist() const
	{
		return _currSolution.dist;
	}

	const std::vector<int>& getCurrSolutionPath() const
	{
		return _currSolution.path;
	}
protected:
	const std::vector<std::vector<int>>* _adjMat = nullptr;
	TSPSolution _currSolution;
	std::mt19937 _gen;
};