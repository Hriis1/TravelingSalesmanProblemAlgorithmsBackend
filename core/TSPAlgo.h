#pragma once
#include <random>

class TSPAlgo
{
public:
	void reseed(unsigned int seed = std::random_device{}())
	{
		_gen.seed(seed);
	}
private:
	std::mt19937 _gen;
};