#pragma once
#include <utility>
#include <string>
#include <array>
#include <vector>
#include <regex>

#include "httplib.h"
#include "json.hpp"

#include "TspApiUtils.h"
#include "TspApiError.h"
#include "../core/TSPUtils.h"
#include "../core/TSPLibParser.h"

#include "../core/TSPGeneticAlgo.h"
#include "../core/TSPMMAS.h"
#include "../core/TSPLKH.h"

namespace tspapiroutes {
	void solveTSP(const httplib::Request& req, httplib::Response& res);

	void getTspInstanceCoords(const httplib::Request& req, httplib::Response& res);

	void getTspCustomFileCoords(const httplib::Request& req, httplib::Response& res);
}