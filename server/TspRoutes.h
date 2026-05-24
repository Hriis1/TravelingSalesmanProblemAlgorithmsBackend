#pragma once
#include <utility>
#include <string>
#include <vector>

#include "httplib.h"
#include "json.hpp"

#include "TspApiUtils.h"
#include "../core/TSPUtils.h"
#include "../core/TSPLibParser.h"

#include "../core/TSPGeneticAlgo.h"
#include "../core/TSPMMAS.h"
#include "../core/TSPLKH.h"


void solveTSP(const httplib::Request& req, httplib::Response& res);