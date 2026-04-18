#pragma once
#include <random>
#include <vector>

#include "TSPSolution.h"

class TSPAlgo
{
public:
	TSPAlgo(unsigned int seed)
		: _gen(seed)
	{
	}

	void reseed(unsigned int seed = std::random_device{}())
	{
		_gen.seed(seed);
	}

	virtual void solve(const std::vector<std::vector<int>>& adjMat) = 0;

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
	TSPSolution _currSolution;
	std::mt19937 _gen;
};