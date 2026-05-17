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

    //Validate the initial tour
    assert(validateTourInternal(startCity));

    //Build segment metadata for the two-level tour representation.
    rebuildTourSegments(startCity);
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

bool TSPLKH::isEdgeInTour(int a, int b) const
{
    return _tourInternal[a].prev == b || _tourInternal[a].next == b;
}

bool TSPLKH::isBetweenInTour(int a, int b, int c) const
{
    const int n = (int)_tourInternal.size();
    assert(a >= 0 && a < n);
    assert(b >= 0 && b < n);
    assert(c >= 0 && c < n);
    assert((int)_citySegment.size() == n);
    assert((int)_cityOffsetInSegment.size() == n);

    //Strict between: endpoints themselves are not considered between.
    if (a == b || b == c || a == c)
        return false;

    auto tourPosition = [&](int city) -> std::pair<int, int>
    {
        const int segmentIndex = _citySegment[city];
        const int offset = _cityOffsetInSegment[city];

        assert(segmentIndex >= 0 && segmentIndex < (int)_tourSegments.size());
        assert(offset >= 0 && offset < (int)_tourSegments[segmentIndex].cities.size());

        const auto& segment = _tourSegments[segmentIndex];

        //If a segment is lazily reversed later, logical order is opposite.
        const int logicalOffset =
            segment.reversed ? (int)segment.cities.size() - 1 - offset : offset;

        return { segment.rank, logicalOffset };
    };

    const auto posA = tourPosition(a);
    const auto posB = tourPosition(b);
    const auto posC = tourPosition(c);

    //If a comes before c in linear tour order, b must be inside that interval.
    if (posA < posC)
        return posA < posB && posB < posC;

    //Otherwise the interval wraps around the end of the tour.
    return posA < posB || posB < posC;
}

long long TSPLKH::calculateInternalTourCost(const std::vector<std::vector<int>>& adjMat, int startCity) const
{
    long long dist = 0;
    int currCity = startCity;

    do
    {
        dist += adjMat[currCity][_tourInternal[currCity].next];
        currCity = _tourInternal[currCity].next;
    } while (currCity != startCity);

    return dist;
}

std::vector<int> TSPLKH::buildOutputPath(int startCity) const
{
    int n = _tourInternal.size();
    assert(n > 0);

    std::vector<int> path(n);

    //add start city
    int currCity = startCity;
    path[0] = startCity;

    //add the other cities
    for (size_t i = 1; i < n; i++)
    {
        currCity = _tourInternal[currCity].next;
        path[i] = currCity;
    }

    return path;
}

void TSPLKH::refreshTourRanks(int startCity)
{
    int currCity = startCity;
    int rank = 0;

    do
    {
        _tourInternal[currCity].rank = rank++;
        currCity = _tourInternal[currCity].next;
    } while (currCity != startCity);
}

void TSPLKH::rebuildTourSegments(int startCity)
{
    const int n = (int)_tourInternal.size();
    assert(n >= 3);
    assert(validateTourInternal(startCity));

    //Use about sqrt(n) cities per segment unless the config overrides it.
    const int targetSegmentSize =
        _config.tourSegmentSize > 0 ? _config.tourSegmentSize : std::max(1, (int)std::sqrt((double)n));

    //Rebuild all segment data from the current linked tour.
    _tourSegments.clear();
    _citySegment.assign(n, -1);
    _cityOffsetInSegment.assign(n, -1);

    int currCity = startCity;
    int visited = 0;

    //Walk the tour once and split consecutive cities into fixed-size blocks.
    while (visited < n)
    {
        LKHTourSegment segment;
        segment.rank = (int)_tourSegments.size();
        segment.reversed = false;
        segment.cities.reserve(targetSegmentSize);

        //Fill one segment with consecutive cities in current tour order.
        while (visited < n && (int)segment.cities.size() < targetSegmentSize)
        {
            const int offset = (int)segment.cities.size();
            segment.cities.push_back(currCity);

            //Remember where each city lives so later order queries are O(1).
            _citySegment[currCity] = segment.rank;
            _cityOffsetInSegment[currCity] = offset;

            currCity = _tourInternal[currCity].next;
            visited++;
        }

        _tourSegments.push_back(segment);
    }

    const int segmentCount = (int)_tourSegments.size();

    //Link the segments themselves into a circular list.
    for (int i = 0; i < segmentCount; i++)
    {
        _tourSegments[i].prev = (i == 0 ? segmentCount - 1 : i - 1);
        _tourSegments[i].next = (i + 1 == segmentCount ? 0 : i + 1);
    }

    assert(validateTourSegments(startCity));
}

bool TSPLKH::validateTourSegments(int startCity) const
{
    const int n = (int)_tourInternal.size();
    if (n < 3)
        return false;

    if (startCity < 0 || startCity >= n)
        return false;

    if ((int)_citySegment.size() != n || (int)_cityOffsetInSegment.size() != n)
        return false;

    if (_tourSegments.empty())
        return false;

    //Each city must appear in exactly one segment.
    std::vector<char> seenCity(n, 0);
    int totalCities = 0;

    for (int segmentIndex = 0; segmentIndex < (int)_tourSegments.size(); segmentIndex++)
    {
        const auto& segment = _tourSegments[segmentIndex];

        //Segment rank should match its index in the segment vector.
        if (segment.rank != segmentIndex)
            return false;

        //Segment links must form a valid circular doubly-linked list.
        if (segment.prev < 0 || segment.prev >= (int)_tourSegments.size())
            return false;

        if (segment.next < 0 || segment.next >= (int)_tourSegments.size())
            return false;

        if (_tourSegments[segment.next].prev != segmentIndex)
            return false;

        if (_tourSegments[segment.prev].next != segmentIndex)
            return false;

        for (int offset = 0; offset < (int)segment.cities.size(); offset++)
        {
            const int city = segment.cities[offset];
            if (city < 0 || city >= n)
                return false;

            //No city may appear in two different segments.
            if (seenCity[city])
                return false;

            seenCity[city] = 1;

            //The city-to-segment lookup must point back to this exact location.
            if (_citySegment[city] != segmentIndex)
                return false;

            if (_cityOffsetInSegment[city] != offset)
                return false;

            totalCities++;
        }
    }

    if (totalCities != n)
        return false;

    for (int city = 0; city < n; city++)
    {
        if (!seenCity[city])
            return false;
    }

    return true;
}

void TSPLKH::applyTwoOptMove(int t1, int t2, int t3, int t4)
{
    assert(_tourInternal[t1].next == t2);
    assert(_tourInternal[t3].prev == t4);

    //Reverse the segment between t2 and t4.
    int currCity = t2;
    while (true)
    {
        std::swap(_tourInternal[currCity].prev, _tourInternal[currCity].next);
        if (currCity == t4)
            break;

        currCity = _tourInternal[currCity].prev;
    }

    //Reconnect the two broken tour edges in their improved arrangement.
    _tourInternal[t1].next = t4;
    _tourInternal[t4].prev = t1;
    _tourInternal[t2].next = t3;
    _tourInternal[t3].prev = t2;

    refreshTourRanks(0);
    rebuildTourSegments(0);

    assert(validateTourInternal(0));
}
