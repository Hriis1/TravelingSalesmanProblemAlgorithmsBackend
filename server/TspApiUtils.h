#pragma once

#include <string>
#include <vector>
#include <utility>
#include <initializer_list>

#include "httplib.h"
#include "json.hpp"

namespace tspapiutils
{
    template <typename Container>
    std::string buildJsonFromPairs(const Container& data)
    {
        nlohmann::json j;

        for (const auto& param : data)
            j[param.first] = param.second;

        return j.dump();
    }

    inline std::string buildJson(const std::vector<std::pair<std::string, nlohmann::json>>& data)
    {
        return buildJsonFromPairs(data);
    }

    inline std::string buildJson(std::initializer_list<std::pair<std::string, nlohmann::json>> data)
    {
        return buildJsonFromPairs(data);
    }

    inline void sendResponse(httplib::Response& res, int status, const std::string& body)
    {
        res.status = status;
        res.set_content(body, "application/json");
    }
}