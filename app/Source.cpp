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
	//testGeneticInstance("ch130.tsp", 500, 100, 0.50f, 0.05f, 100, true, 10);

	//Run MMAS
	//testMMASRand(100, 1000, 20, 500, 2, 3, 0.1, 100, false);
	//testMMASInstance("brg180.tsp", 500, 2, 3, 0.1, 100, 10);

	//Run LKH
	testLKHInstance("fl1400.tsp", 1);


	std::cin.get();
	return 0;
}