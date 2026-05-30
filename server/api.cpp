#include <iostream>

#include "TspRoutes.h"

int main()
{
    httplib::Server server;

    //Get instance coords
    server.Get("/getTspInstanceCoords", tspapiroutes::getTspInstanceCoords);

    //Solve tsp
    server.Post("/solveTSP", tspapiroutes::solveTSP);

    std::cout << "Server running on http://localhost:8080\n";
    server.listen("0.0.0.0", 8080);
}