#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <random>
#include <algorithm>
#include <climits>

namespace TSPUtils
{

    enum class TspDatasetType {
        RandomUniform,
        Circular,
        ClusteredDeceptive
    };

    int distRoundedInt(double x1, double y1, double x2, double y2)
    {
        const double dx = x1 - x2;
        const double dy = y1 - y2;
        return (int)std::lround(std::sqrt(dx * dx + dy * dy));
    }

    std::vector<std::vector<int>> generateTspAdjMatrix(unsigned int numCities, TspDatasetType type, int planeSize = 1000, unsigned int seed = std::random_device{}())
    {
        if (numCities < 2) return {};

        std::mt19937 rng(seed);
        std::uniform_real_distribution<double> uni(0.0, (double)planeSize);

        std::vector<std::pair<double, double>> pts;
        pts.reserve(numCities);

        if (type == TspDatasetType::RandomUniform)
        {
            for (unsigned int i = 0; i < numCities; i++)
                pts.push_back({ uni(rng), uni(rng) });
        }
        else if (type == TspDatasetType::Circular)
        {
            // Points on a circle + small noise
            const double cx = planeSize * 0.5;
            const double cy = planeSize * 0.5;
            const double R = planeSize * 0.40;

            std::normal_distribution<double> noise(0.0, planeSize * 0.01); // 1% noise

            for (unsigned int i = 0; i < numCities; i++)
            {
                const double angle = (2.0 * M_PI * (double)i) / (double)numCities;
                double x = cx + R * std::cos(angle) + noise(rng);
                double y = cy + R * std::sin(angle) + noise(rng);
                x = std::max(0.0, std::min(x, (double)planeSize)); //clamp
                y = std::max(0.0, std::min(y, (double)planeSize)); //clamp
                pts.push_back({ x, y });
            }

            // Optional: shuffle so the “correct” tour isn't trivially 0..n-1
            std::shuffle(pts.begin(), pts.end(), rng);
        }
        else // ClusteredDeceptive
        {
            // Make K cluster centers far apart, then sample points near centers.
            // This tends to create hard "cluster ordering" decisions.
            const unsigned int K = std::max(2u, (unsigned int)std::sqrt((double)numCities / 2.0));
            std::vector<std::pair<double, double>> centers;
            centers.reserve(K);

            // Place centers in a rough grid to ensure separation
            const unsigned int grid = (unsigned int)std::ceil(std::sqrt((double)K));
            const double step = planeSize / (double)(grid + 1);

            unsigned int placed = 0;
            for (unsigned int gy = 1; gy <= grid && placed < K; gy++)
            {
                for (unsigned int gx = 1; gx <= grid && placed < K; gx++)
                {
                    // jitter center a bit
                    std::normal_distribution<double> jitter(0.0, step * 0.08);
                    double x = gx * step + jitter(rng);
                    double y = gy * step + jitter(rng);
                    x = std::max(0.0, std::min(x, (double)planeSize)); //clamp
                    y = std::max(0.0, std::min(y, (double)planeSize)); //clamp
                    centers.push_back({ x, y });
                    placed++;
                }
            }

            // Tight cluster radius
            const double sigma = planeSize * 0.02; // 2% of plane
            std::normal_distribution<double> clNoise(0.0, sigma);

            for (unsigned int i = 0; i < numCities; i++)
            {
                const auto& c = centers[i % K];
                double x = c.first + clNoise(rng);
                double y = c.second + clNoise(rng);
                x = std::max(0.0, std::min(x, (double)planeSize)); //clamp
                y = std::max(0.0, std::min(y, (double)planeSize)); //clamp
                pts.push_back({ x, y });
            }

            std::shuffle(pts.begin(), pts.end(), rng);
        }

        // Build adjacency matrix
        std::vector<std::vector<int>> adj(numCities, std::vector<int>(numCities, 0));
        for (unsigned int i = 0; i < numCities; i++)
        {
            for (unsigned int j = i + 1; j < numCities; j++)
            {
                int d = distRoundedInt(pts[i].first, pts[i].second, pts[j].first, pts[j].second);
                adj[i][j] = adj[j][i] = d;
            }
        }
        return adj;
    }

    int bruteForceOptimal(const std::vector<std::vector<int>>& adjMatrix, int startCity = 0)
    {
        const int n = (int)adjMatrix.size();

        // Edge cases
        if (n == 0) 
            return 0;

        if (n == 1)
            return 0;

        if (startCity < 0 || startCity >= n) 
            return -1;

        for (const auto& row : adjMatrix)
            if ((int)row.size() != n)
                return -1;

        // Build list of cities excluding startCity
        std::vector<int> perm;
        perm.reserve(n - 1);
        for (int i = 0; i < n; ++i)
            if (i != startCity)
                perm.push_back(i);

        // perm is already sorted by construction
        long long best = std::numeric_limits<long long>::max();

        do
        {
            long long dist = 0;

            int prev = startCity;
            for (int city : perm)
            {
                dist += adjMatrix[prev][city];
                prev = city;

                // small prune (optional)
                if (dist >= best) break;
            }

            if (dist < best)
            {
                dist += adjMatrix[prev][startCity]; // close the tour

                if (dist < best) 
                    best = dist;
            }
        } while (std::next_permutation(perm.begin(), perm.end()));

        if (best > std::numeric_limits<int>::max()) 
            return INT_MAX; // avoid overflow on return

        return (int)best;
    }

    int nearestNeighborDistance(const std::vector<std::vector<int>>& adjMatrix, int startCity = 0)
    {
        const int n = adjMatrix.size();

        // Edge cases
        if (n == 0)
            return 0;

        if (n == 1)
            return 0;

        if (startCity < 0 || startCity >= n)
            return -1;

        for (const auto& row : adjMatrix)
            if (row.size() != n)
                return -1;

        std::vector<char> visited(n, 0);
        visited[startCity] = 1;

        int currCity = startCity;
        int totalDist = 0;

        for (int step = 1; step < n; ++step)
        {
            int nextCity = -1;
            int bestDist = INT_MAX;

            for (int city = 0; city < n; ++city)
            {
                if (visited[city])
                    continue;

                const int d = adjMatrix[currCity][city];
                if (d < bestDist)
                {
                    bestDist = d;
                    nextCity = city;
                }
            }

            if (nextCity < 0)
                return -1;

            visited[nextCity] = 1;
            totalDist += bestDist;
            currCity = nextCity;
        }

        totalDist += adjMatrix[currCity][startCity];
        return totalDist;
    }

    std::vector<int> nearestNeighborPath(const std::vector<std::vector<int>>& adjMatrix, int startCity = 0)
    {
        const int n = adjMatrix.size();

        // Edge cases
        if (n == 0) 
            return {};

        if (n == 1)
        {
            if (startCity != 0) 
                return {}; // invalid start for 1-city matrix

            return { 0 }; //just start city
        }

        // Validate start city
        if (startCity < 0 || startCity >= n) 
            return {};

        // Validate matrix is n x n
        for (const auto& row : adjMatrix)
            if (row.size() != n) 
                return {};

        std::vector<char> visited(n, 0);
        visited[startCity] = 1;

        std::vector<int> path;
        path.reserve(n);
        path.push_back(startCity);

        int currCity = startCity;

        for (int step = 1; step < n; ++step)
        {
            int nextCity = -1;
            int bestDist = INT_MAX;

            for (int city = 0; city < n; ++city)
            {
                if (visited[city]) 
                    continue;

                const int d = adjMatrix[currCity][city];
                if (d < bestDist)
                {
                    bestDist = d;
                    nextCity = city;
                }
            }

            if (nextCity < 0) 
                return {}; // shouldn't happen

            visited[nextCity] = 1;
            path.push_back(nextCity);
            currCity = nextCity;
        }

        return path;
    }

}
