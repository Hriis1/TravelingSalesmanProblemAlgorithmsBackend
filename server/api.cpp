#include <iostream>

#include "TspRoutes.h"

int main()
{
    httplib::Server server;

    server.Post("/solveTSP", tspApiRoutes::solveTSP);

    std::cout << "Server running on http://localhost:8080\n";
    server.listen("0.0.0.0", 8080);
}