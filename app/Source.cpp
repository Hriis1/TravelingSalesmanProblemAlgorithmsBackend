#include <iostream>
#include <vector>
#include <iomanip>

#include "../core/TSPGeneticAlgo.h"
#include "../core/TSPUtils.h"
#include "../core/TestTSPAlgos.h"
#include "../core/TSPLibParser.h"

int main() 
{
	//Run genetic algo
	//testGeneticAlgoRand(100, 1000, 10, 500, 100, 0.90f, 0.05f, 100, false, true);

	//Run MMAS
	testMMASRand(500, 1000, 2, 500, 2, 3, 0.1, 100, false);

	std::cin.get();
	return 0;
}