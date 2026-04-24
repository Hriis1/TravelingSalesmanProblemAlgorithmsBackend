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
	void twoOpt(std::vector<int>& path, const std::vector<std::vector<int>>& adjMat)
	{
		int n = path.size();
		bool improved = true;

		while (improved) {
			improved = false;

			for (int i = 1; i < n - 1; i++) {
				for (int j = i + 1; j < n; j++) {

					int a = path[i - 1];
					int b = path[i];
					int c = path[j];
					int d = path[(j + 1) % n];

					int oldDist = adjMat[a][b] + adjMat[c][d];
					int newDist = adjMat[a][c] + adjMat[b][d];

					if (newDist < oldDist) {
						//reverse segment [i, j]
						std::reverse(path.begin() + i, path.begin() + j + 1);
						improved = true;
					}
				}
			}
		}
	}

	int calculatePathDist(const std::vector<int>& path, const std::vector<std::vector<int>>& adjMat)
	{
		int n = path.size();
		if (n == 0) return 0;

		int dist = 0;

		for (int i = 0; i < n - 1; i++) {
			dist += adjMat[path[i]][path[i + 1]];
		}

		// return to start
		dist += adjMat[path[n - 1]][path[0]];

		return dist;
	}
protected:
	TSPSolution _currSolution;
	std::mt19937 _gen;
};