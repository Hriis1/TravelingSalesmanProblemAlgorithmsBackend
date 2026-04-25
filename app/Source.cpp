#include <iostream>
#include <vector>
#include <iomanip>

#include "../core/TSPGeneticAlgo.h"
#include "../core/TSPUtils.h"
#include "../core/TestTSPAlgos.h"

int main() 
{
	//Run genetic algo
	//testGeneticAlgo(100, 1000, 10, 500, 100, 0.90f, 0.05f, 100, false, true);

	//Run MMAS
	testMMAS(100, 1000, 10, 1000, 1, 2, 0.5, 100, false);

	std::cin.get();
	return 0;
}