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
    assert(n >= 3);
    assert(startCity >= 0 && startCity < n);

    std::vector<char> seenCity(n, 0);
    std::vector<char> seenRank(n, 0);

    int curr = startCity;
    for (int step = 0; step < n; step++)
    {
        assert(curr >= 0 && curr < n);
        assert(!seenCity[curr]);
        seenCity[curr] = 1;

        const auto& node = _tourInternal[curr];
        assert(node.prev >= 0 && node.prev < n);
        assert(node.next >= 0 && node.next < n);
        assert(node.rank >= 0 && node.rank < n);
        assert(!seenRank[node.rank]);
        seenRank[node.rank] = 1;

        assert(_tourInternal[node.next].prev == curr);
        assert(_tourInternal[node.prev].next == curr);

        curr = node.next;
    }

    assert(curr == startCity);

    for (int i = 0; i < n; i++)
    {
        assert(seenCity[i]);
        assert(seenRank[i]);
    }

    return true;
}
