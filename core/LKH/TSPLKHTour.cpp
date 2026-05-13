#include "../TSPLKH.h"

void TSPLKH::buildInitialTourNN(const std::vector<std::vector<int>>& adjMat, int startCity)
{
    int n = adjMat.size();

    //Add start city
    std::vector<char> visited(n, 0);
    visited[startCity] = 1;
    addCityToTour(startCity, -1, 0);
    int currCity = startCity;

    //Choose next city
    for (int step = 1; step < n; ++step)
    {
        int nextCity = -1;
        int bestDist = INT_MAX;

        for (int city = 0; city < n; ++city)
        {
            if (visited[city])
                continue;

            const int d = adjMat[currCity][city];
            if (d < bestDist)
            {
                bestDist = d;
                nextCity = city;
            }
        }

        //Add next city and hook it to curr city
        visited[nextCity] = 1;
        addCityToTour(nextCity, currCity, step);
        _tourInternal[currCity].next = nextCity;
        currCity = nextCity;
    }

    //Finish the cycle
    _tourInternal[currCity].next = startCity;
    _tourInternal[startCity].prev = currCity;
}

void TSPLKH::addCityToTour(int city, int prev, int rank)
{
    _tourInternal[city].prev = prev;
    _tourInternal[city].rank = rank;
}
