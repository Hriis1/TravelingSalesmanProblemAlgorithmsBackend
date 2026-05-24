#include "TspRoutes.h"

//empty namespace for internal functions
namespace
{

}

void solveTSP(const httplib::Request& req, httplib::Response& res)
{
    nlohmann::json bodyJson;

    try {
        bodyJson = nlohmann::json::parse(req.body);
    }
    catch (...) { //Body is invalid json
        tspapiutils::sendResponse(res, 400,
            tspapiutils::buildJson({ {"success", false}, {"error", "Invalid JSON body"} }));
        return;
    }

    //Invalid TSP instance param
    if (!bodyJson.contains("instance") || !bodyJson["instance"].is_string()) {
        tspapiutils::sendResponse(res, 400,
            tspapiutils::buildJson({ {"success", false}, {"error", "Missing or invalid field: instance"} }));
        return;
    }

    //Invalid Algorithm param 
    if (!bodyJson.contains("algorithm") || !bodyJson["algorithm"].is_string()) {
        tspapiutils::sendResponse(res, 400,
            tspapiutils::buildJson({ {"success", false}, {"error", "Missing or invalid field: algorithm"} }));
        return;
    }

    //Bulld final response
    std::string instance = bodyJson["instance"];
    std::string algorithm = bodyJson["algorithm"];
    auto response = tspapiutils::buildJson({
        {"success", true},
        {"instance", instance},
        {"algorithm", algorithm}
        });

    //Send final response
    tspapiutils::sendResponse(res, 200, response);
}