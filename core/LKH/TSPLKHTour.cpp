#include "../TSPLKH.h"
#include <algorithm>

void TSPLKH::buildInitialTourNN(const std::vector<std::vector<int>>& adjMat, int startCity)
{
    int n = adjMat.size();
    std::vector<char> visited(n, 0);
    std::vector<int> nnPath(n);

    //Add start city
    int currCity = startCity;
    visited[currCity] = 1;
    nnPath[0] = currCity;

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
        nnPath[step] = nextCity;
        currCity = nextCity;
    }

    //Build segment metadata for the two-level tour representation.
    rebuildTourSegmentsFromPath(nnPath);

    rebuildTourInternalFromSegments(); // temporary cache until _tourInternal is fully phased out
}

void TSPLKH::addCityToTour(int city, int prev, int rank)
{
    _tourInternal[city].prev = prev;
    _tourInternal[city].rank = rank;
}

bool TSPLKH::validateTourInternal(int startCity) const
{
    const int n = (int)_citySegment.size();
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
    return getPrevInTour(a) == b || getNextInTour(a) == b;
}

int TSPLKH::getPrevInTour(int city) const
{
    const int n = (int)_citySegment.size();
    assert(city >= 0 && city < n);
    assert((int)_citySegment.size() == n);
    assert((int)_cityOffsetInSegment.size() == n);

    const int segmentIndex = _citySegment[city];
    const int offset = _cityOffsetInSegment[city];

    assert(segmentIndex >= 0 && segmentIndex < (int)_tourSegments.size());
    const auto& segment = _tourSegments[segmentIndex];
    assert(offset >= 0 && offset < (int)segment.cities.size());

    //Inside a normal segment, previous city is one physical slot left.
    if (!segment.reversed && offset > 0)
        return segment.cities[offset - 1];

    //Inside a reversed segment, logical previous city is one physical slot right.
    if (segment.reversed && offset + 1 < (int)segment.cities.size())
        return segment.cities[offset + 1];

    //At a segment boundary, previous city is the logical last city of the previous segment.
    const auto& prevSegment = _tourSegments[segment.prev];
    return prevSegment.reversed ? prevSegment.cities.front() : prevSegment.cities.back();
}

int TSPLKH::getNextInTour(int city) const
{
    const int n = (int)_citySegment.size();
    assert(city >= 0 && city < n);
    assert((int)_citySegment.size() == n);
    assert((int)_cityOffsetInSegment.size() == n);

    const int segmentIndex = _citySegment[city];
    const int offset = _cityOffsetInSegment[city];

    assert(segmentIndex >= 0 && segmentIndex < (int)_tourSegments.size());
    const auto& segment = _tourSegments[segmentIndex];
    assert(offset >= 0 && offset < (int)segment.cities.size());

    //Inside a normal segment, next city is one physical slot right.
    if (!segment.reversed && offset + 1 < (int)segment.cities.size())
        return segment.cities[offset + 1];

    //Inside a reversed segment, logical next city is one physical slot left.
    if (segment.reversed && offset > 0)
        return segment.cities[offset - 1];

    //At a segment boundary, next city is the logical first city of the next segment.
    const auto& nextSegment = _tourSegments[segment.next];
    return nextSegment.reversed ? nextSegment.cities.back() : nextSegment.cities.front();
}

bool TSPLKH::isBetweenInTour(int a, int b, int c) const
{
    const int n = (int)_citySegment.size();
    assert(a >= 0 && a < n);
    assert(b >= 0 && b < n);
    assert(c >= 0 && c < n);
    assert((int)_citySegment.size() == n);
    assert((int)_cityOffsetInSegment.size() == n);

    //Strict between: endpoints themselves are not considered between.
    if (a == b || b == c || a == c)
        return false;

    //Convert a city into a comparable logical tour position.
    //The first value is the segment order, the second is the position inside that segment.
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
        const int nextCity = getNextInTour(currCity);
        dist += adjMat[currCity][nextCity];
        currCity = nextCity;
    } while (currCity != startCity);

    return dist;
}

std::vector<int> TSPLKH::buildOutputPath(int startCity) const
{
    int n = (int)_citySegment.size();
    assert(n > 0);

    std::vector<int> path(n);

    //add start city
    int currCity = startCity;
    path[0] = startCity;

    //add the other cities
    for (size_t i = 1; i < n; i++)
    {
        currCity = getNextInTour(currCity);
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
        currCity = getNextInTour(currCity);
    } while (currCity != startCity);
}

void TSPLKH::rebuildTourSegments(int startCity)
{
    const int n = (int)_citySegment.size();
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
    const int n = (int)_citySegment.size();
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

bool TSPLKH::validateTourStructure(int startCity) const
{
    const int n = (int)_citySegment.size();
    if (n < 3)
        return false;

    if (startCity < 0 || startCity >= n)
        return false;

    //First validate the raw segment containers and city lookup tables.
    if (!validateTourSegments(startCity))
        return false;

    std::vector<char> seenCity(n, 0);
    int currCity = startCity;

    //Walk the tour using the segment-aware next/prev helpers.
    for (int step = 0; step < n; step++)
    {
        if (currCity < 0 || currCity >= n)
            return false;

        if (seenCity[currCity])
            return false;

        seenCity[currCity] = 1;

        const int nextCity = getNextInTour(currCity);
        const int prevCity = getPrevInTour(currCity);

        //Next and previous links must agree with each other.
        if (getPrevInTour(nextCity) != currCity)
            return false;

        if (getNextInTour(prevCity) != currCity)
            return false;

        currCity = nextCity;
    }

    //After visiting n unique cities, the next step must close the cycle.
    if (currCity != startCity)
        return false;

    for (int city = 0; city < n; city++)
    {
        if (!seenCity[city])
            return false;
    }

    return true;
}

void TSPLKH::rebuildTourSegmentsFromPath(const std::vector<int>& path)
{
    const int n = (int)path.size();
    assert(n >= 3);

    //Use about sqrt(n) cities per segment unless the config overrides it.
    const int targetSegmentSize =
        _config.tourSegmentSize > 0 ? _config.tourSegmentSize : std::max(1, (int)std::sqrt((double)n));

    std::vector<char> seenCity(n, 0);

    _tourSegments.clear();

    for (int pathIndex = 0; pathIndex < n;)
    {
        LKHTourSegment segment;
        segment.rank = (int)_tourSegments.size();
        segment.reversed = false;
        segment.cities.reserve(targetSegmentSize);

        //Copy one consecutive block from the path into this segment.
        while (pathIndex < n && (int)segment.cities.size() < targetSegmentSize)
        {
            const int city = path[pathIndex++];
            assert(city >= 0 && city < n);
            assert(!seenCity[city]);

            seenCity[city] = 1;
            segment.cities.push_back(city);
        }

        _tourSegments.push_back(segment);
    }

    //Rebuild segment links and city -> segment lookup arrays.
    rebuildSegmentIndexes();

    assert(validateTourSegments(path[0]));
}

void TSPLKH::rebuildSegmentIndexes()
{
    const int n = (int)_citySegment.size();
    const int segmentCount = (int)_tourSegments.size();

    _citySegment.assign(n, -1);
    _cityOffsetInSegment.assign(n, -1);

    for (int segmentIndex = 0; segmentIndex < segmentCount; segmentIndex++)
    {
        auto& segment = _tourSegments[segmentIndex];

        //Keep the vector order, segment rank, and segment links in sync.
        segment.rank = segmentIndex;
        segment.prev = segmentIndex == 0 ? segmentCount - 1 : segmentIndex - 1;
        segment.next = segmentIndex + 1 == segmentCount ? 0 : segmentIndex + 1;

        //City offsets are physical positions inside segment.cities.
        for (int offset = 0; offset < (int)segment.cities.size(); offset++)
        {
            const int city = segment.cities[offset];
            _citySegment[city] = segmentIndex;
            _cityOffsetInSegment[city] = offset;
        }
    }
}

std::vector<int> TSPLKH::getSegmentCitiesInTourOrder(const LKHTourSegment& segment) const
{
    //Return the segment contents in logical tour order.
    //If the segment is marked reversed, physical storage is read backwards.
    std::vector<int> cities = segment.cities;
    if (segment.reversed)
        std::reverse(cities.begin(), cities.end());

    return cities;
}

void TSPLKH::splitSegmentBeforeCity(int city)
{
    assert(city >= 0 && city < (int)_citySegment.size());

    const int segmentIndex = _citySegment[city];
    auto logicalCities = getSegmentCitiesInTourOrder(_tourSegments[segmentIndex]);

    const auto it = std::find(logicalCities.begin(), logicalCities.end(), city);
    assert(it != logicalCities.end());

    const int splitOffset = (int)(it - logicalCities.begin());
    if (splitOffset == 0)
        return;

    //Split the old segment into the cities before 'city' and the cities starting at 'city'.
    LKHTourSegment before;
    before.cities.assign(logicalCities.begin(), logicalCities.begin() + splitOffset);

    LKHTourSegment after;
    after.cities.assign(logicalCities.begin() + splitOffset, logicalCities.end());

    //The split pieces are stored directly in logical order, so their lazy flag is reset.
    _tourSegments[segmentIndex] = before;
    _tourSegments.insert(_tourSegments.begin() + segmentIndex + 1, after);
    rebuildSegmentIndexes();
}

void TSPLKH::splitSegmentAfterCity(int city)
{
    assert(city >= 0 && city < (int)_citySegment.size());

    const int segmentIndex = _citySegment[city];
    auto logicalCities = getSegmentCitiesInTourOrder(_tourSegments[segmentIndex]);

    const auto it = std::find(logicalCities.begin(), logicalCities.end(), city);
    assert(it != logicalCities.end());

    const int splitOffset = (int)(it - logicalCities.begin()) + 1;
    if (splitOffset == (int)logicalCities.size())
        return;

    //Split the old segment into the cities through 'city' and the cities after it.
    LKHTourSegment before;
    before.cities.assign(logicalCities.begin(), logicalCities.begin() + splitOffset);

    LKHTourSegment after;
    after.cities.assign(logicalCities.begin() + splitOffset, logicalCities.end());

    //The split pieces are stored directly in logical order, so their lazy flag is reset.
    _tourSegments[segmentIndex] = before;
    _tourSegments.insert(_tourSegments.begin() + segmentIndex + 1, after);
    rebuildSegmentIndexes();
}

void TSPLKH::rotateTourSegmentsToFront(int segmentIndex)
{
    assert(segmentIndex >= 0 && segmentIndex < (int)_tourSegments.size());

    //A circular tour has no fixed first segment, so rotating keeps the same tour.
    std::rotate(_tourSegments.begin(), _tourSegments.begin() + segmentIndex, _tourSegments.end());
    rebuildSegmentIndexes();
}

void TSPLKH::rebuildTourInternalFromSegments()
{
    const int n = (int)_citySegment.size();
    std::vector<int> tourOrder;
    tourOrder.reserve(n);

    //Materialize the logical segment order back into one city order.
    //This keeps older code that reads _tourInternal.prev/next correct.
    for (const auto& segment : _tourSegments)
    {
        const auto cities = getSegmentCitiesInTourOrder(segment);
        tourOrder.insert(tourOrder.end(), cities.begin(), cities.end());
    }

    assert((int)tourOrder.size() == n);

    for (int rank = 0; rank < n; rank++)
    {
        const int city = tourOrder[rank];
        _tourInternal[city].prev = tourOrder[rank == 0 ? n - 1 : rank - 1];
        _tourInternal[city].next = tourOrder[rank + 1 == n ? 0 : rank + 1];
        _tourInternal[city].rank = rank;
    }
}

void TSPLKH::flipTourPath(int firstCity, int lastCity)
{
    const int n = (int)_citySegment.size();
    assert(firstCity >= 0 && firstCity < n);
    assert(lastCity >= 0 && lastCity < n);
    assert(validateTourInternal(0));

    if (firstCity == lastCity)
        return;

    //Make the flip boundaries align with segment boundaries.
    splitSegmentBeforeCity(firstCity);
    splitSegmentAfterCity(lastCity);

    //After splitting, firstCity is the first city of its segment.
    const int firstSegment = _citySegment[firstCity];
    rotateTourSegmentsToFront(firstSegment);

    //After rotation, the path from firstCity forward starts at segment 0.
    const int lastSegment = _citySegment[lastCity];

    //Reverse the segment order for the flipped path.
    std::reverse(_tourSegments.begin(), _tourSegments.begin() + lastSegment + 1);

    //Each segment inside the flipped path also has its internal direction reversed.
    for (int i = 0; i <= lastSegment; i++)
        _tourSegments[i].reversed = !_tourSegments[i].reversed;

    //Only segment metadata is updated here. Materialization is done by the caller if needed.
    rebuildSegmentIndexes();

    assert(validateTourStructure(0));
}

void TSPLKH::applyTwoOptMove(int t1, int t2, int t3, int t4)
{
    assert(getNextInTour(t1) == t2);
    assert(getPrevInTour(t3) == t4);

    //Reverse the tour path that remains between the two deleted edges.
    //Because the tour is cyclic, this also creates the new boundary edges:
    //t1 -> t4 and t2 -> t3.
    flipTourPath(t2, t4);

    //Materialize once so older code that still reads _tourInternal remains correct.
    rebuildTourInternalFromSegments();

    assert(validateTourStructure(0));
    assert(validateTourInternal(0));
}
