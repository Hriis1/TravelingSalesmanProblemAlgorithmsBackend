#include <iostream>
#include <utility>
#include <string>
#include <vector>

#include "httplib.h"
#include "json.hpp"

template <typename Container>
std::string buildJsonFromPairs(const Container& data)
{
    nlohmann::json j;

    for (const auto& param : data)
        j[param.first] = param.second;

    return j.dump();
}

std::string buildJson(const std::vector<std::pair<std::string, nlohmann::json>>& data)
{
    return buildJsonFromPairs(data);
}

std::string buildJson(std::initializer_list<std::pair<std::string, nlohmann::json>> data)
{
    return buildJsonFromPairs(data);
}

void sendResponse(httplib::Response& res, int status, const std::string& body)
{
    res.status = status;
    res.set_content(body, "application/json");
}

void solveTSP(const httplib::Request& req, httplib::Response& res)
{
    nlohmann::json bodyJson;

    try {
        bodyJson = nlohmann::json::parse(req.body);
    }
    catch (...) { //Body is invalid json
        sendResponse(res, 400,
            buildJson({ {"success", false}, {"error", "Invalid JSON body"} }));
        return;
    }

    //Invalid TSP instance param
    if (!bodyJson.contains("instance") || !bodyJson["instance"].is_string()) {
        sendResponse(res, 400,
            buildJson({ {"success", false}, {"error", "Missing or invalid field: instance"} }));
        return;
    }

    //Invalid Algorithm param 
    if (!bodyJson.contains("algorithm") || !bodyJson["algorithm"].is_string()) {
        sendResponse(res, 400,
            buildJson({ {"success", false}, {"error", "Missing or invalid field: algorithm"} }));
        return;
    }

    //Bulld final response
    std::string instance = bodyJson["instance"];
    std::string algorithm = bodyJson["algorithm"];
    auto response = buildJson({
        {"success", true},
        {"instance", instance},
        {"algorithm", algorithm}
        });

    //Send final response
    sendResponse(res, 200, response);
}

int main()
{
    httplib::Server server;

    server.Post("/solveTSP", solveTSP);

    std::cout << "Server running on http://localhost:8080\n";
    server.listen("0.0.0.0", 8080);
}