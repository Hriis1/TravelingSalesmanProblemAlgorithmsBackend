#include "httplib.h"
#include <iostream>

int main()
{
    httplib::Server server;

    server.Get("/test", [](const httplib::Request& req, httplib::Response& res) {
        res.set_content("{\"message\":\"API works\"}", "application/json");
        });

    std::cout << "Server running on http://localhost:8080\n";
    server.listen("0.0.0.0", 8080);
}