#pragma once

#include <iostream>
#include <vector>
#include <iomanip>
#include <chrono>
#include <string>

#include "../core/TSPGeneticAlgo.h"
#include "../core/TSPMMAS.h"
#include "../core/TSPUtils.h"
#include "../core/TSPLibParser.h"

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

int testTSPAlgoRand(int nCities, int planeSize, int nRuns, bool doBruteForce, TSPAlgo* tspSolver, const std::string& algoName)
{
	//Input adj matrix
	std::vector<std::vector<int>> adjMat;

	//avgs
	int genAvg = 0;
	int nearestNeighborAvg = 0;
	int optimalAvg = 0;
	long long msPassedTotal = 0;

	//Generate the matrix
	std::cout << "Running " << algoName  << " and nearest neighbor " << nRuns << " times..." << std::endl;
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
		auto start = std::chrono::high_resolution_clock::now();
		tspSolver->solve(adjMat);
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		long long msPassed = duration.count();
		msPassedTotal += msPassed;

		//Output time passed
		std::cout << "TSP Solved in: " << msPassed << " ms" << std::endl;

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

		std::cout << std::endl << std::endl;
	}

	//Do avg and display
	genAvg /= nRuns;
	nearestNeighborAvg /= nRuns;
	float pDecrease = (float)genAvg / nearestNeighborAvg;
	std::cout << std::endl << std::endl << algoName << " avg: " << genAvg << std::endl << "Nearest neighbor avg: " << nearestNeighborAvg << std::endl;
	std::cout << "On avg " << algoName << " tour is " << (pDecrease * 100) << "% of the NN tour length" << std::endl;
	std::cout << "Avg time to solve TSP: " << msPassedTotal / nRuns << " ms" << std::endl;


	if (doBruteForce) {
		optimalAvg /= nRuns;
		pDecrease = (float)genAvg / optimalAvg;
		std::cout << std::endl << std::endl << "Optimal avg: " << optimalAvg << std::endl;
		std::cout << "% decreese in tour lenght from OPTIMAL: " << (pDecrease * 100) << "%";
	}

	return 1;
}

int testTSPAlgoInstance(const std::string& tspInstance, TSPAlgo* tspSolver, const std::string& algoName)
{
	const int nRuns = 10;
	try
	{
		//Init the instance
		auto instance = TSPLibParser::parseFile("repo/tsplib-master/" + tspInstance);

		if (instance.optimalDist == -1)
		{
			std::cout << "No known optimal distance for instance: " << instance.name << std::endl;
			return 0;
		}

		if (instance.adjMat.size() == 0 || instance.adjMat.size() != instance.adjMat[0].size())
		{
			std::cout << "Invalid TSPLIB instance matrix!" << std::endl;
			return 0;
		}

		long long distTotal = 0;
		long long msPassedTotal = 0;
		double gapTotal = 0.0;
		int nearestNeighborDist = TSPUtils::nearestNeighborDistance(instance.adjMat, 0, false);
		double gapNN = 100.0 * (double)(nearestNeighborDist - instance.optimalDist) / (double)instance.optimalDist;


		std::cout << "Running TSP instance " << algoName << " for " << instance.name << " " << nRuns << " times..." << std::endl;
		std::cout << "Dimension: " << instance.dimension << std::endl;
		std::cout << "Optimal dist: " << instance.optimalDist << std::endl << std::endl;

		for (int i = 0; i < nRuns; i++)
		{


			auto start = std::chrono::high_resolution_clock::now();
			tspSolver->solve(instance.adjMat);
			auto end = std::chrono::high_resolution_clock::now();

			auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
			long long msPassed = duration.count();
			int dist = tspSolver->getCurrSolutionDist();
			double gap = 100.0 * (double)(dist - instance.optimalDist) / (double)instance.optimalDist;

			distTotal += dist;
			msPassedTotal += msPassed;
			gapTotal += gap;

			std::cout << "Run " << (i + 1) << ":" << std::endl;
			std::cout << "Dist: " << dist << std::endl;
			std::cout << "Time: " << msPassed << " ms" << std::endl;
			std::cout << "Gap from optimal: " << std::fixed << std::setprecision(2) << gap << "%" << std::endl;
			std::cout << std::endl;
		}

		std::cout << std::endl << "Averages:" << std::endl;
		std::cout << "Avg dist: " << (double)distTotal / nRuns << std::endl;
		std::cout << "Nearest neighbor dist: " << nearestNeighborDist << std::endl;
		std::cout << "Avg time: " << (double)msPassedTotal / nRuns << " ms" << std::endl;
		std::cout << "Optimal dist: " << instance.optimalDist << std::endl << std::endl;
		std::cout << "Avg gap from optimal: " << std::fixed << std::setprecision(2) << gapTotal / nRuns << "%" << std::endl;
		std::cout << "Nearest neighbor gap from optimal: " << std::fixed << std::setprecision(2) << gapNN << "%" << std::endl;
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		return 0;
	}
	

	return 1;
}

int testGeneticAlgoRand(int nCities, int planeSize, int nRuns, int ng, int npop, float pc, float pm, int nnoimpr, bool doBruteForce, bool initWithNN)
{
	//Init solver
	TSPGeneticAlgo tsp = TSPGeneticAlgo(ng, npop, nnoimpr, pc, pm, initWithNN);

	//Solve and display results
	return testTSPAlgoRand(nCities, planeSize, nRuns, doBruteForce, &tsp, "Genetic algorithm");
}

int testGeneticInstance(const std::string& instance, int ng, int npop, float pc, float pm, int nnoimpr, bool initWithNN)
{
	//Init solver
	TSPGeneticAlgo tsp = TSPGeneticAlgo(ng, npop, nnoimpr, pc, pm, initWithNN);

	//Solve and display results
	return testTSPAlgoInstance(instance, &tsp, "Genetic algorithm");
}

int testMMASRand(int nCities, int planeSize, int nRuns, int nIters, double alpha, double beta, double rho, int nnoimpr, bool doBruteForce)
{
	//Init solver
	TSPMMAS tsp = TSPMMAS(nIters, alpha, beta, rho, nnoimpr);

	//Solve and display results
	return testTSPAlgoRand(nCities, planeSize, nRuns, doBruteForce, &tsp, "MMAS");
}

int testMMASInstance(const std::string& instance, int nIters, double alpha, double beta, double rho, int nnoimpr) 
{
	//Init solver
	TSPMMAS tsp = TSPMMAS(nIters, alpha, beta, rho, nnoimpr);

	//Solve and display results
	return testTSPAlgoInstance(instance, &tsp, "MMAS");
}

