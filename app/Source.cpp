#include <iostream>
#include <vector>
#include <iomanip>

#include "../core/TSPGeneticAlgo.h"
#include "../core/TSPUtils.h"
#include "../core/TestTSPAlgos.h"

int main() 
{
	//Run genetic algo
	//testGeneticAlgoRand(100, 1000, 10, 500, 100, 0.90f, 0.05f, 100, false, true);

	//Run MMAS
	//testMMASRand(100, 1000, 20, 500, 2, 3, 0.1, 100, false);
	testMMASInstance("a280.tsp", 500, 2, 3, 0.1, 100);

	std::cin.get();
	return 0;
}