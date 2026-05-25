#pragma once
#include <string>
#include <stdexcept>

class TspApiError : public std::runtime_error
{
public:
    int statusCode;

    TspApiError(int code, const std::string& message)
        : std::runtime_error(message), statusCode(code) {}
};