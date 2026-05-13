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

    assert(validateTourInternal(startCity));
}

void TSPLKH::addCityToTour(int city, int prev, int rank)
{
    _tourInternal[city].prev = prev;
    _tourInternal[city].rank = rank;
}

bool TSPLKH::validateTourInternal(int startCity) const
{
    const int n = (int)_tourInternal.size();
    if (n < 3)
        return false;

    if (startCity < 0 || startCity >= n)
        return false;

    std::vector<char> seenCity(n, 0);
    std::vector<char> seenRank(n, 0);

    int curr = startCity;
    for (int step = 0; step < n; step++)
    {
        if (curr < 0 || curr >= n)
            return false;

        if (seenCity[curr])
            return false;

        seenCity[curr] = 1;

        const auto& node = _tourInternal[curr];
        if (node.prev < 0 || node.prev >= n)
            return false;

        if (node.next < 0 || node.next >= n)
            return false;

        if (node.rank < 0 || node.rank >= n)
            return false;

        if (seenRank[node.rank])
            return false;

        seenRank[node.rank] = 1;

        if (_tourInternal[node.next].prev != curr)
            return false;

        if (_tourInternal[node.prev].next != curr)
            return false;

        curr = node.next;
    }

    if (curr != startCity)
        return false;

    for (int i = 0; i < n; i++)
    {
        if (!seenCity[i])
            return false;

        if (!seenRank[i])
            return false;
    }

    return true;
}
