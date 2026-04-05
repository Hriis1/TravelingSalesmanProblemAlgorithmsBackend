#pragma once

#include <vector>
#include <algorithm>
#include <random>
#include <cassert>
#include <climits>

struct Genome
{
public:

	//vars
	int dist = INT_MAX;
	std::vector<int> path;

	//funcs
	Genome(int numCities) 
	{
		path.reserve(numCities);
	}

    void calculateDist(const std::vector<std::vector<int>>& adjMat)
    {
        //Reset dist
        dist = 0;

        //Sz of path
        int N = path.size();

        if (N == 0)
            return;

        //Asert path size == size of adj matrix
        assert(N == adjMat.size());

        //Calc dist
        for (size_t i = 0; i < path.size() - 1; i++)
            dist += adjMat[path[i]][path[i + 1]];

        // add return to start
        dist += adjMat[path[N - 1]][path[0]];
    }

    void generate(const std::vector<std::vector<int>>& adjMat, std::mt19937& gen)
    {
        int N = adjMat.size();

        path.clear();

        //Empty adj matrix edge case
        if (N == 0)
        {
            dist = 0;
            return;
        }

        // 1. Fill with all cities
        for (int i = 0; i < N; ++i)
            path.emplace_back(i);

        // 2. Shuffle everything except first element since city 0 is the start
        std::shuffle(path.begin() + 1, path.end(), gen);

        // 3. Compute distance
        dist = 0;

        for (int i = 0; i < N - 1; ++i)
            dist += adjMat[path[i]][path[i + 1]];

        // add return to start
        dist += adjMat[path[N - 1]][path[0]];
    }
};