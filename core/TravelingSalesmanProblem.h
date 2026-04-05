#pragma once
#include <climits>
#include <array>
#include <cassert>

#include "Genome.h"
#include "TSPUtils.h"


class TravelingSalesmanProblem
{
public:
	TravelingSalesmanProblem(const std::vector<std::vector<int>>& adjMat, int ng, int npop, int nnoimpr, float pc, float pm, unsigned int seed = std::random_device{}())
		:_adjMat(&adjMat), _NG(ng), _NPOP(npop), _NNOIMPR(nnoimpr), _PC(pc), _PM(pm), _currSolution(adjMat.size()), _gen(seed)
	{
		_currSolution.dist = INT_MAX;

		assert(_NPOP >= 3);
	}

	void reinit(const std::vector<std::vector<int>>& adjMat, int ng, int nnoimpr, int npop, float pc, float pm)
	{
		_adjMat = &adjMat;
		_NG = ng;
		_NPOP = npop;
		_NNOIMPR = nnoimpr;
		_PC = pc;
		_PM = pm;
		_currSolution = Genome(adjMat.size());
		_currSolution.dist = INT_MAX;

		assert(_NPOP >= 3);
	}

	void reseed(unsigned int seed = std::random_device{}())
	{
		_gen.seed(seed);
	}

	void solve(bool useNNIn1stGen = false)
	{
		//Init variables
		int noImproveCounter = _NNOIMPR;
		int currBestDist = INT_MAX;
		int nParticipants = std::min(4, _NPOP);
		std::vector<Genome> currGen(_NPOP, Genome(_adjMat->size()));
		std::vector<Genome> nextGen;
		nextGen.reserve(_NPOP);
		std::uniform_real_distribution<float> probDist(0.0f, 1.0f);

		//Generate the 1st generation
		size_t i = 0;
		if (useNNIn1stGen)
		{
			//Use nearest neighbor as a child in the initial generation
			currGen[i].path = TSPUtils::nearestNeighborPath(*_adjMat, 0);
			currGen[i].calculateDist(*_adjMat);
			i++;
		}
		//Generare the rest randomly
		for (; i < _NPOP; i++)
		{
			currGen[i].generate(*_adjMat, _gen);
		}

		//Get the idxes of the 2 best solutions
		int bestIdx1 = -1;
		int bestIdx2 = -1;

		//helper func for finding the idxes of the 2 best solutions
		auto findBestTwo = [&](const std::vector<Genome>& generation)
		{
			bestIdx1 = -1;
			bestIdx2 = -1;

			for (size_t i = 0; i < static_cast<int>(generation.size()); ++i)
			{
				if (bestIdx1 == -1 || generation[i].dist < generation[bestIdx1].dist)
				{
					bestIdx2 = bestIdx1;
					bestIdx1 = i;
				}
				else if (i != bestIdx1 &&
					(bestIdx2 == -1 || generation[i].dist < generation[bestIdx2].dist))
				{
					bestIdx2 = i;
				}
			}
		};

		//find the idxes of the 2 best solutions of this gen
		findBestTwo(currGen);
		currBestDist = currGen[bestIdx1].dist;

		//Generate new _NG - 1 new generations or stop when there hasn't been an improvement for __NNOIMPR generations
		for (size_t i = 1; i < _NG; i++)
		{
			//Reset next gen
			nextGen.clear();
			nextGen.reserve(_NPOP);

			//Elitism - copy the top 2 genomes of curr gen
			//best
			nextGen.emplace_back(currGen[bestIdx1]);
			//2nd best
			nextGen.emplace_back(currGen[bestIdx2]);

			//Fill remaining elements using tournament selection, crossover and mutation
			for (size_t n = 2; n < _NPOP; n++)
			{
				//select parent 1 idx using tournament selection
				int p1Idx = tournamentSelection(currGen, nParticipants);

				//based on the _PC probabilty copy p1 or crossover with p2
				if (probDist(_gen) < _PC)
				{
					// do crossover
					int p2Idx = tournamentSelection(currGen, nParticipants);
					auto childPath = edgeRecombinationCrossover(currGen[p1Idx].path, currGen[p2Idx].path);
					nextGen.emplace_back(_adjMat->size());
					nextGen[n].path = std::move(childPath);
				}
				else
				{
					// just copy p1
					nextGen.emplace_back(currGen[p1Idx]);
				}

				//base on the _PM probability mutate the new child
				if (probDist(_gen) < _PM)
				{
					displacementMutation(nextGen[n].path);
				}

				//Calculate the dist of the new child
				nextGen[n].calculateDist(*_adjMat);
			}

			//Find the idxes of the 2 best solutions of this gen
			findBestTwo(nextGen);

			//A new best solution was found
			if (currBestDist > nextGen[bestIdx1].dist) 
			{
				noImproveCounter = _NNOIMPR;
				currBestDist = nextGen[bestIdx1].dist;
			}
			else //A new best solution was not found
			{
				//End algo if counter becomes 0
				if (--noImproveCounter <= 0)
				{
					_currSolution =  nextGen[bestIdx1];
					return;
				}
			}

			//Move curr gen to next gen
			currGen = std::move(nextGen);
		}

		//Set the best solution
		_currSolution = currGen[bestIdx1];
	}

	//getters
	int getCurrSolutionDist() const
	{
		return _currSolution.dist;
	}

	const std::vector<int>& getCurrSolutionPath() const
	{
		return _currSolution.path;
	}

private:

	std::vector<int> edgeRecombinationCrossover(const std::vector<int>& p1, const std::vector<int>& p2) 
	{

		//Asert that both parents have correct size for correctness
		assert(p2.size() == p1.size());

		//Size of parent/child
		int N = p1.size();

		//child vec
		std::vector<int> child(N);

		//Empty parent
		if (N == 0)
			return child;

		//Adjacency matrix and num neighbors for all vertecies
		std::vector<std::array<int, 4>> adjMatrix(N);
		std::vector<uint8_t> deg(N, 0);

		//Build the adj matrix from both parents
		// helper: add neighbor if not already present
		auto addUnique = [&](int city, int neigh)
		{
			if (neigh == city) return;

			for (uint8_t t = 0; t < deg[city]; ++t)
				if (adjMatrix[city][t] == neigh)
					return; // duplicate

			// assumes deg[city] < 4
			adjMatrix[city][deg[city]++] = neigh;
		};

		//p1, dont check for duplicates
		for (size_t i = 0; i < N; i++)
		{
			int city = p1[i];

			//Get the left/right neighbors
			int leftCity = p1[(i - 1 + N) % N];
			int rightCity = p1[(i + 1) % N];

			// Add to adjacency list and increment deg
			adjMatrix[city][deg[city]++] = leftCity;
			adjMatrix[city][deg[city]++] = rightCity;
		}

		// p2, check for duplicates
		for (size_t i = 0; i < N; ++i)
		{
			int city = p2[i];

			//Get the left/right neighbors
			int leftCity = p2[(i - 1 + N) % N];
			int rightCity = p2[(i + 1) % N];

			// Add to adjacency list and increment deg, checking for duplicates
			addUnique(city, leftCity);
			addUnique(city, rightCity);

			//Asert deg < 4 for correctness
			assert(deg[city] <= 4);
		}

		// Build the child from adj matrix
		// helper: remove "rem" from neighbor list of "city"
		auto removeNeighbor = [&](int city, int rem)
		{
			for (uint8_t i = 0; i < deg[city]; ++i)
			{
				if (adjMatrix[city][i] == rem)
				{
					// swap with last and shrink
					adjMatrix[city][i] = adjMatrix[city][deg[city] - 1];
					deg[city]--;
					return;
				}
			}
		};

		// helper: after choosing city "chosen", remove it from all its neighbors' lists
		auto removeChosenFromAll = [&](int chosen)
		{
			// chosen currently has up to deg[chosen] neighbors
			uint8_t d = deg[chosen];
			for (uint8_t i = 0; i < d; ++i)
			{
				int u = adjMatrix[chosen][i];
				removeNeighbor(u, chosen);
			}
			//deg[chosen] = 0; // clear chosen's list
		};

		//used cities
		std::vector<char> used(N, 0);

		// --- O(1) fallback support: keep a dynamic list of unused cities ---
		std::vector<int> remaining;
		remaining.reserve(N);
		std::vector<int> posInRemaining(N, -1);

		for (size_t c = 0; c < N; ++c)
		{
			posInRemaining[c] = (int)remaining.size();
			remaining.push_back(c);
		}

		// helper: remove a city from the remaining list in O(1)
		auto removeFromRemaining = [&](int city)
		{
			int idx = posInRemaining[city];
			if (idx < 0) return; // already removed

			int last = remaining.back();
			remaining[idx] = last;
			posInRemaining[last] = idx;

			remaining.pop_back();
			posInRemaining[city] = -1;
		};

		// first city
		child[0] = p1[0];
		used[child[0]] = 1;
		removeFromRemaining(child[0]);
		removeChosenFromAll(child[0]);

		// ERX construction loop
		for (size_t i = 1; i < N; ++i)
		{
			int lastCity = child[i - 1];

			// collect best candidates among available (not used) neighbors of lastCity
			int bestDeg = INT_MAX;
			std::vector<int> candidates;
			candidates.reserve(4);

			for (uint8_t n = 0; n < deg[lastCity]; ++n)
			{
				int neigh = adjMatrix[lastCity][n];
				if (used[neigh]) continue;

				int d = (int)deg[neigh];
				if (d < bestDeg)
				{
					bestDeg = d;
					candidates.clear();
					candidates.push_back(neigh);
				}
				else if (d == bestDeg)
				{
					candidates.push_back(neigh);
				}
			}

			int nextCity = -1;

			if (!candidates.empty())
			{
				// tie-break randomly among min-degree candidates
				std::uniform_int_distribution<int> pick(0, (int)candidates.size() - 1);
				nextCity = candidates[pick(_gen)];
			}
			else
			{
				// fallback: pick a random unused city (O(1) now)
				std::uniform_int_distribution<int> pick(0, (int)remaining.size() - 1);
				nextCity = remaining[pick(_gen)];
			}

			child[i] = nextCity;
			used[nextCity] = 1;
			removeFromRemaining(nextCity);

			// critical ERX step: remove chosen city from all adjacency lists
			removeChosenFromAll(nextCity);
		}

		return child;
	}

	void displacementMutation(std::vector<int>& child)
	{
		int N = child.size();

		if (N < 3) 
			return; // with index 0 fixed, there's nothing meaningful to displace

		// Pick i, j with 1 <= i <= j <= N-1
		std::uniform_int_distribution<int> distI(1, N - 1);
		int i = distI(_gen);

		std::uniform_int_distribution<int> distJ(i, N - 1);
		int j = distJ(_gen);

		// If the segment is the entire mutable tail, there's nowhere to insert it
		const int segLen = j - i + 1;
		if (segLen == N - 1)
			return;

		// Pick k in [1..N-1] such that k is NOT in [i..j] and k != i-1 (to avoid no-op)
		std::uniform_int_distribution<int> distK(1, N - 1);
		int k = 0;
		while(true)
		{
			k = distK(_gen);

			if (k >= i && k <= j) 
				continue;     // must be outside the segment

			break;
		}

		// Perform displacement using O(N) rotate
		// insert at index k
		auto b = child.begin();

		if (k < i)
		{
			// ... k | (k+1..i-1) | [i..j] | ...
			// -> ... k | [i..j] | (k+1..i-1) | ...
			std::rotate(b + k, b + i, b + (j + 1));
		}
		else // k > j
		{
			// ... [i..j] | (j+1..k) | ...
			// -> ... (j+1..k) | [i..j] | ...
			std::rotate(b + i, b + (j + 1), b + k);
		}
	}

	int tournamentSelection(const std::vector<Genome>& generation, int numParticipants)
	{
		//generation should have as many elements as defined in _NPOP
		assert(_NPOP == generation.size());

		//Cant have more participants than num genomes in a generation
		assert(_NPOP >= numParticipants);

		//Generate random numParticipants participants and return the best of them
		int currBestIdx = std::uniform_int_distribution<int>(0, _NPOP - 1)(_gen);
		int currBestDist = generation[currBestIdx].dist;
		for (size_t i = 1; i < numParticipants; i++)
		{
			int rIdx = std::uniform_int_distribution<int>(0, _NPOP - 1)(_gen);
			int currDist = generation[rIdx].dist;
			if (currDist < currBestDist) 
			{
				currBestDist = currDist;
				currBestIdx = rIdx;
			}
		}

		return currBestIdx;
	}

	const std::vector<std::vector<int>>* _adjMat = nullptr;
	int _NG = 0;
	int _NPOP = 0;
	int _NNOIMPR = 0;
	float _PC = 0.0f;
	float _PM = 0.0f;
	Genome _currSolution;
	std::mt19937 _gen;
};
