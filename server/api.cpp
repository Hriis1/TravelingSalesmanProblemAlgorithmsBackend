#include <iostream>
#include <utility>
#include <string>
#include <vector>

#include "httplib.h"
#include "json.hpp"

std::string buildJson(const std::vector<std::pair<std::string, nlohmann::json>>& data)
{
    nlohmann::json j;

    for (const auto& param : data)
        j[param.first] = param.second;

    return j.dump();
}

void solveTSP(const httplib::Request& req, httplib::Response& res)
{
    auto body = buildJson({
        {"success", true},
        {"message", "API works"},
        {"distance", 10466}
    });

    res.status = 200;
    res.set_content(body, "application/json");
}

int main()
{
    httplib::Server server;

    server.Get("/solveTSP", solveTSP);

    std::cout << "Server running on http://localhost:8080\n";
    server.listen("0.0.0.0", 8080);
}