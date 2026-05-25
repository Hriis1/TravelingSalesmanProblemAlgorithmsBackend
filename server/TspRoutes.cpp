#include "TspRoutes.h"

//empty namespace for internal functions
namespace
{
    void solveTSPInstance(const httplib::Request& req, httplib::Response& res, const nlohmann::json& reqBody, TSPAlgo& tspSolver)
    {
        //Init instance
        std::string sourceDir = TSPUtils::getSourceDir();
        const std::string& tspInstance = reqBody["instance"];
        std::string path = sourceDir + "/../tsplib-master/" + tspInstance + ".tsp";
        auto instance = TSPLibParser::parseFile(path);

        //Instance errors
        if (instance.adjMat.size() == 0 || instance.adjMat.size() != instance.adjMat[0].size())
        {
            throw std::runtime_error("Invalid TSPLIB instance matrix!");
        }

        //Run the solver
        int nearestNeighborDist = TSPUtils::nearestNeighborDistance(instance.adjMat, 0);
        tspSolver.solve(instance.adjMat);

        //Send response
        tspapiutils::sendResponse(res, 200,
            tspapiutils::buildJson({ {"success", true}, {"path", tspSolver.getCurrSolutionPath()}, {"nCities", instance.adjMat.size()},
                {"dist", tspSolver.getCurrSolutionDist()}, {"nnDist", nearestNeighborDist}, {"optimalDist", instance.optimalDist} }));
    }
}

void solveTSP(const httplib::Request& req, httplib::Response& res)
{
    nlohmann::json bodyJson;

    try
    {
        try {
            bodyJson = nlohmann::json::parse(req.body);
        }
        catch (...) { //Body is invalid json
            throw std::runtime_error("Invalid JSON body");
        }

        //Invalid Algorithm param 
        if (!bodyJson.contains("algorithm") || !bodyJson["algorithm"].is_string()) {
            throw std::runtime_error("Missing or invalid field: algorithm");
        }

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
            throw std::runtime_error("algorithm not recognized");
        }

        //Solve instance
        if (bodyJson.contains("instance") && bodyJson["instance"].is_string()) {

            solveTSPInstance(req, res, bodyJson, *tspSolver);
            return;
        }
    }
    catch (const std::exception& e)
    {
        //Catch errors
        std::cerr << e.what() << std::endl;
        tspapiutils::sendResponse(res, 400,
            tspapiutils::buildJson({ {"success", false}, {"error", e.what()} }));
        return;
    }
}