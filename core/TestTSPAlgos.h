#pragma once

#include <iostream>
#include <vector>
#include <iomanip>

#include "../core/TSPGeneticAlgo.h"
#include "../core/TSPUtils.h"

void outputPath(const std::vector<int>& path)
{
	const int sz = path.size();

	if (sz == 0)
		return;

	//Output path
	for (size_t i = 0; i < sz; i++)
	{
		std::cout << path[i] << " -> ";
	}

	//Return to start city
	if (sz > 1)
	{
		std::cout << path[0];
	}

	//New line
	std::cout << std::endl;
}

void printMatrix(const std::vector<std::vector<int>>& mat)
{
	for (const auto& row : mat)
	{
		for (const auto& val : row)
		{
			std::cout << std::setw(3) << val << " ";
		}
		std::cout << std::endl;
	}
}

int testGeneticAlgo()
{
	//If brute force should be done
	bool doBruteForce = false;

	//If initial gen should be inited with NN
	bool initWithNN = true;

	//Number of runs to do
	const int nRuns = 10;

	//Input adj matrix
	std::vector<std::vector<int>> adjMat;

	//avgs
	int genAvg = 0;
	int nearestNeighborAvg = 0;
	int optimalAvg = 0;

	//Generate the matrix
	std::cout << "Running genetic algorithm and nearest neighbor " << nRuns << " times..." << std::endl;
	for (size_t i = 0; i < nRuns; i++)
	{
		std::cout << "Genereting matrix..." << std::endl;
		adjMat = TSPUtils::generateTspAdjMatrix(100, TSPUtils::TspDatasetType::RandomUniform);
		std::cout << "Matrix generated!" << std::endl;

		//Print the matrix
		//std::cout << "Matrix:" << std::endl;
		//printMatrix(adjMat);

		//if matrix is not square => invalid input
		if (adjMat.size() == 0 || (adjMat.size() != adjMat[0].size()))
		{
			std::cout << "Invalid input!" << std::endl;
			return -1;
		}

		//Nubmer of cities
		const int nCities = adjMat.size();

		//Set up params;
		int ng = 0; //number of generation
		int npop = 0; //number of children in a generation
		int nnoimpr = 0; //after how many generations with no improvement to stop
		float pc = 0.90f; //probability of crossover
		float pm = 0.05f; //probability of mutation

		//set up params based on num cities
		ng = 500;
		npop = 100;
		nnoimpr = 100;

		//Solve - unseeded
		TSPGeneticAlgo tsp = TSPGeneticAlgo(ng, npop, nnoimpr, pc, pm, initWithNN);
		tsp.solve(adjMat);

		//Output path and dist
		int genDist = tsp.getCurrSolutionDist();
		int nearestNeighborDist = TSPUtils::nearestNeighborDistance(adjMat, 0);
		std::cout << "Path: " << std::endl;
		outputPath(tsp.getCurrSolutionPath());
		std::cout << "Dist: " << genDist << std::endl;
		std::cout << "Nearest-neighbor dist (start city 0): " << nearestNeighborDist << std::endl;

		//Accumulate avg
		genAvg += genDist;
		nearestNeighborAvg += nearestNeighborDist;

		if (doBruteForce) {
			int optimalDist = TSPUtils::bruteForceOptimal(adjMat, 0);
			std::cout << "Optimal dist (start city 0): " << optimalDist << std::endl;
			optimalAvg += optimalDist;
		}
	}

	//Do avg and display
	genAvg /= nRuns;
	nearestNeighborAvg /= nRuns;
	float pDecrease = (float)genAvg / nearestNeighborAvg;
	std::cout << std::endl << std::endl << "Genetic algo avg: " << genAvg << std::endl << "Nearest neighbor avg: " << nearestNeighborAvg << std::endl;
	std::cout << "% decreese in tour lenght from NN: " << (pDecrease * 100) << "%" << std::endl;

	if (doBruteForce) {
		optimalAvg /= nRuns;
		pDecrease = (float)genAvg / optimalAvg;
		std::cout << std::endl << std::endl << "Optimal avg: " << optimalAvg << std::endl;
		std::cout << "% decreese in tour lenght from OPTIMAL: " << (pDecrease * 100) << "%";
	}

	return 0;
}