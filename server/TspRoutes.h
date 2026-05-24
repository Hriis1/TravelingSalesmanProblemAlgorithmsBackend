#pragma once
#include <utility>
#include <string>
#include <vector>

#include "TspApiUtils.h"
#include "httplib.h"
#include "json.hpp"

void solveTSP(const httplib::Request& req, httplib::Response& res);