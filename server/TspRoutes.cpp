#include "TspRoutes.h"

//empty namespace for internal functions
namespace
{
    std::vector<std::vector<double>> jsonCitiesToCoords(const nlohmann::json& cities)
    {
        int n = static_cast<int>(cities.size());

        std::vector<std::vector<double>> coords;
        coords.reserve(n);

        for (const auto& city : cities)
        {
            if (!city.is_array() || city.size() != 2 ||
                !city[0].is_number() || !city[1].is_number())
            {
                throw TspApiError(400, "Each city must be an array: [x, y]");
            }

            coords.push_back({
                city[0].get<double>(),
                city[1].get<double>()
                });
        }

        return coords;
    }

    void solveTSPInstance(const httplib::Request& req, httplib::Response& res, const nlohmann::json& reqBody, TSPAlgo& tspSolver)
    {
        //Init instance
        std::string sourceDir = TSPUtils::getSourceDir();
        const std::string& tspInstance = reqBody["instance"];
        std::string path = sourceDir + "/../tsplib-master/" + tspInstance + ".tsp";
        auto instance = [&]() {
            try
            {
                return TSPLibParser::parseFile(path);
            }
            catch (const std::exception& e)
            {
                throw TspApiError(400, e.what());
            }
        }();

        //Instance errors
        if (instance.adjMat.size() == 0 || instance.adjMat.size() != instance.adjMat[0].size())
            throw TspApiError(400, "Invalid TSPLIB instance matrix!");

        //Run the solver
        int nearestNeighborDist = TSPUtils::nearestNeighborDistance(instance.adjMat, 0);
        tspSolver.solve(instance.adjMat);

        //Send response
        tspapiutils::sendResponse(res, 200,
            tspapiutils::buildJson({ {"success", true}, {"path", tspSolver.getCurrSolutionPath()}, {"nCities", instance.adjMat.size()},
                {"dist", tspSolver.getCurrSolutionDist()}, {"nnDist", nearestNeighborDist}, {"optimalDist", instance.optimalDist} }));
    }

    void solveTSPCustom(const httplib::Request& req, httplib::Response& res, const nlohmann::json& reqBody, TSPAlgo& tspSolver)
    {
        const int minCities = 5;
        const int maxCities = 3000;

        const auto& custom = reqBody["customTSP"];

        //Num cities exists
        if (!custom.contains("numCities"))
            throw TspApiError(400, "Missing field: customTSP.numCities");

        //cities exist
        if (!custom.contains("cities"))
            throw TspApiError(400, "Missing field: customTSP.cities");

        const auto& numCitiesJson = reqBody["customTSP"]["numCities"];
        const auto& cities = reqBody["customTSP"]["cities"];

        //numCities is an int
        if (!numCitiesJson.is_number_integer())
            throw TspApiError(400, "numCities must be an integer");

        //numCities between min and max cities
        int numCities = numCitiesJson.get<int>();
        if (numCities < minCities || numCities > maxCities)
            throw TspApiError(400, "numCities must be between 5 and 3000");

        //count of cities must match numCities
        if (!cities.is_array() || cities.size() != numCities)
            throw TspApiError(400, "cities must be an array of the same size as numCities");

        //parse cities to adj matrix
        std::vector<std::vector<double>> coords = jsonCitiesToCoords(cities);
        std::vector<std::vector<int>> adjMat = TSPUtils::buildAdjMatrixFromCoords(coords);

        //Run the solver
        int nearestNeighborDist = TSPUtils::nearestNeighborDistance(adjMat, 0);
        tspSolver.solve(adjMat);

        //Send response
        tspapiutils::sendResponse(res, 200,
            tspapiutils::buildJson({ {"success", true}, {"path", tspSolver.getCurrSolutionPath()}, {"nCities", adjMat.size()},
                {"dist", tspSolver.getCurrSolutionDist()}, {"nnDist", nearestNeighborDist}, {"optimalDist", -1} }));
    }
}

void solveTSP(const httplib::Request& req, httplib::Response& res)
{
    try
    {
        nlohmann::json bodyJson;

        try {
            bodyJson = nlohmann::json::parse(req.body);
        }
        catch (...) { //Body is invalid json
            throw TspApiError(400, "Invalid JSON body");
        }

        //Invalid Algorithm param 
        if (!bodyJson.contains("algorithm") || !bodyJson["algorithm"].is_string())
            throw TspApiError(400, "Missing or invalid field: algorithm");

        //Init solver
        std::unique_ptr<TSPAlgo> tspSolver;
        if (bodyJson["algorithm"] == "genetic")
        {
            int ng = 500;
            int npop = 100;
            int nnoimpr = 100;
            float pc = 0.90f;
            float pm = 0.1f;
            tspSolver = std::make_unique<TSPGeneticAlgo>(ng, npop, nnoimpr, pc, pm, true);
        }
        else if (bodyJson["algorithm"] == "mmas")
        {
            int nIters = 500;
            double alpha = 2;
            double beta = 3;
            double rho = 0.1;
            int nnoimpr = 200;
            tspSolver = std::make_unique<TSPMMAS>(nIters, alpha, beta, rho, nnoimpr);
        }
        else if (bodyJson["algorithm"] == "lkh")
        {
            LKHConfig config;
            tspSolver = std::make_unique<TSPLKH>(std::move(config));
        }
        else
        {
            throw TspApiError(400, "Algorithm not recognized");
        }

        
        //Check for instance or custom
        if (bodyJson.contains("instance") && bodyJson["instance"].is_string()) { //is instance
            //Solve instance
            solveTSPInstance(req, res, bodyJson, *tspSolver);
            return;
        }
        else if (bodyJson.contains("customTSP") && bodyJson["customTSP"].is_object()) //is custom
        {
            //Solve custom tsp
            solveTSPCustom(req, res, bodyJson, *tspSolver);
            return;
        }
        else //No instance or custom sent in req
        {
            throw TspApiError(400, "No valid TSP was sent for solving");
        }
    }
    catch (const TspApiError& e)
    {
        //Handled exceptions
        tspapiutils::sendResponse(res, e.statusCode,
            tspapiutils::buildJson({ {"success", false}, {"error", e.what()} }));
    }
    catch (const std::exception& e)
    {
        //Unhandled exceptions
        tspapiutils::sendResponse(res, 500,
            tspapiutils::buildJson({ {"success", false}, {"error", e.what()} }));
    }
}