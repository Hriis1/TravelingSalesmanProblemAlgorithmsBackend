#pragma once

#include <iostream>
#include <vector>
#include <iomanip>

#include "../core/TSPGeneticAlgo.h"
#include "../core/TSPMMAS.h"
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

int testTSPAlso(int nCities, int planeSize, int nRuns, bool doBruteForce, TSPAlgo* tspSolver)
{
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
		//Generate the matrix
		std::cout << "Genereting matrix..." << std::endl;
		adjMat = TSPUtils::generateTspAdjMatrix(nCities, TSPUtils::TspDatasetType::RandomUniform, planeSize);
		std::cout << "Matrix generated!" << std::endl;

		//Print the matrix
		//std::cout << "Matrix:" << std::endl;
		//printMatrix(adjMat);

		//if matrix is not square => invalid input
		if (adjMat.size() == 0 || (adjMat.size() != adjMat[0].size()))
		{
			std::cout << "Invalid input!" << std::endl;
			return 0;
		}


		//Solve - unseeded
		tspSolver->solve(adjMat);

		//Output path and dist
		int genDist = tspSolver->getCurrSolutionDist();
		int nearestNeighborDist = TSPUtils::nearestNeighborDistance(adjMat, 0);
		std::cout << "Path: " << std::endl;
		outputPath(tspSolver->getCurrSolutionPath());
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

	return 1;
}

int testMMAS(int nCities, int planeSize, int nRuns, int nIters, double alpha, double beta, double rho, int nnoimpr, bool doBruteForce)
{
	//Init solver
	TSPMMAS tsp = TSPMMAS(nCities, nIters, alpha, beta, rho, nnoimpr);

	//Solve and display results
	return testTSPAlso(nCities, planeSize, nRuns, doBruteForce, &tsp);
}

int testGeneticAlgo(int nCities, int planeSize, int nRuns, int ng, int npop, float pc, float pm, int nnoimpr, bool doBruteForce, bool initWithNN)
{
	//Init solver
	TSPGeneticAlgo tsp = TSPGeneticAlgo(ng, npop, nnoimpr, pc, pm, initWithNN);

	//Solve and display results
	return testTSPAlso(nCities, planeSize, nRuns, doBruteForce, &tsp);
}