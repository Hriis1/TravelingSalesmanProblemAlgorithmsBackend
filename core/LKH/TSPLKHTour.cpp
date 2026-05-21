#include "../TSPLKH.h"
#include <algorithm>
#include <numeric>

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
}

void TSPLKH::buildInitialTourWalk(int startCity)
{
    const int n = (int)_candidates.size();
    assert(n >= 3);
    assert(startCity >= 0 && startCity < n);

    std::vector<char> chosen(n, 0);
    std::vector<int> next(n);
    std::vector<int> prev(n);

    //LKH keeps nodes in a circular list and moves each selected successor
    //directly after the current node. This list reproduces that behavior.
    for (int city = 0; city < n; city++)
    {
        next[city] = city + 1 == n ? 0 : city + 1;
        prev[city] = city == 0 ? n - 1 : city - 1;
    }

    std::vector<int> path;
    path.reserve(n);

    int current = startCity;
    chosen[current] = 1;
    path.push_back(current);

    for (int count = 1; count < n; count++)
    {
        std::vector<int> alternatives;

        //Default LKH WALK case D: prefer unchosen candidate edges.
        for (const auto& candidate : _candidates[current])
        {
            const int to = candidate.to;
            if (to >= 0 && to < n && !chosen[to])
                alternatives.push_back(to);
        }

        int selected = -1;
        if (!alternatives.empty())
        {
            //If several candidate successors are possible, LKH chooses one at random.
            std::uniform_int_distribution<int> pick(0, (int)alternatives.size() - 1);
            selected = alternatives[pick(_gen)];
        }
        else
        {
            //Default LKH WALK case E: take the next unchosen node in the current ring.
            selected = next[current];
            while (chosen[selected])
                selected = next[selected];
        }

        //Move the selected node to follow current in the circular list, like LKH Follow().
        if (next[current] != selected)
        {
            next[prev[selected]] = next[selected];
            prev[next[selected]] = prev[selected];

            prev[selected] = current;
            next[selected] = next[current];
            prev[next[current]] = selected;
            next[current] = selected;
        }

        current = selected;
        chosen[current] = 1;
        path.push_back(current);
    }

    rebuildTourSegmentsFromPath(path);
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

void TSPLKH::applyKSwapKick(int k)
{
    const int n = (int)_citySegment.size();
    assert(validateTourStructure(0));

    //LKH's K-swap kick is defined for K >= 4.
    if (k < 4 || n < 4)
        return;

    k = std::min(k, n);

    std::vector<int> path = buildOutputPath(0);

    //Sample K distinct break positions in the current tour.
    std::vector<int> positions(n);
    std::iota(positions.begin(), positions.end(), 0);
    for (int i = 0; i < k; i++)
    {
        std::uniform_int_distribution<int> pick(i, n - 1);
        const int selected = pick(_gen);
        std::swap(positions[i], positions[selected]);
    }

    std::vector<int> breaks(positions.begin(), positions.begin() + k);
    std::sort(breaks.begin(), breaks.end());

    std::vector<int> kickedPath;
    kickedPath.reserve(n);

    auto appendBlock = [&](int blockIndex)
    {
        //Block i starts after break i and ends at the next break city.
        const int start = (breaks[blockIndex] + 1) % n;
        const int end = breaks[(blockIndex + 1) % k];

        int index = start;
        while (true)
        {
            kickedPath.push_back(path[index]);
            if (index == end)
                break;
            index = (index + 1) % n;
        }
    };

    //LKH reconnects s[(i + 2) % K] to old successor of s[i].
    //In block form this means: block 0, then K-1, K-2, ..., 1.
    appendBlock(0);
    for (int blockIndex = k - 1; blockIndex >= 1; blockIndex--)
        appendBlock(blockIndex);

    assert((int)kickedPath.size() == n);

    //Store the kicked tour back in the segmented representation.
    rebuildTourSegmentsFromPath(kickedPath);
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

void TSPLKH::flipTourPath(int firstCity, int lastCity)
{
    const int n = (int)_citySegment.size();
    assert(firstCity >= 0 && firstCity < n);
    assert(lastCity >= 0 && lastCity < n);

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

void TSPLKH::flipPathAndReconnect(int a, int b, int c, int d)
{
    const int n = (int)_citySegment.size();
    assert(a >= 0 && a < n);
    assert(b >= 0 && b < n);
    assert(c >= 0 && c < n);
    assert(d >= 0 && d < n);

    //This primitive represents:
    //remove (a,b) and (d,c), add (a,d) and (b,c).
    assert(getNextInTour(a) == b);
    assert(getNextInTour(d) == c);

    //The transformation is exactly a reversal of the forward path b..d.
    //Because the tour is cyclic, reversing that interval reconnects:
    //a -> d at the left boundary and b -> c at the right boundary.
    flipTourPath(b, d);

    assert(validateTourStructure(0));
}
