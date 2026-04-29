#pragma once

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

struct TSPLibInstance
{
    std::string name;
    std::string type;
    std::string edgeWeightType;
    std::string edgeWeightFormat;
    int dimension = 0;
    int optimalDist = -1;
    std::vector<std::vector<int>> adjMat;
};

class TSPLibParser
{
private:
    struct Point
    {
        double x = 0.0;
        double y = 0.0;
    };

    enum class Section
    {
        None,
        NodeCoord,
        EdgeWeight
    };

public:
    static TSPLibInstance parseFile(const std::string& filePath, int maxDimensionForAdjMatrix = 5000)
    {
        std::ifstream in(filePath);
        if (!in)
            throw std::runtime_error("Could not open TSPLIB file: " + filePath);

        TSPLibInstance instance;
        std::vector<Point> points;
        std::vector<int> edgeWeights;

        Section section = Section::None;
        std::string line;

        while (std::getline(in, line))
        {
            line = trim(line);
            if (line.empty())
                continue;

            const std::string upperLine = toUpper(line);
            if (upperLine == "EOF")
                break;

            if (upperLine == "NODE_COORD_SECTION")
            {
                section = Section::NodeCoord;
                continue;
            }

            if (upperLine == "EDGE_WEIGHT_SECTION")
            {
                section = Section::EdgeWeight;
                continue;
            }

            if (endsWith(upperLine, "_SECTION") || upperLine == "DISPLAY_DATA_SECTION")
            {
                section = Section::None;
                continue;
            }

            if (section == Section::NodeCoord)
            {
                readPointLine(line, instance.dimension, points);
                continue;
            }

            if (section == Section::EdgeWeight)
            {
                readIntLine(line, edgeWeights);
                continue;
            }

            readHeaderLine(line, instance);
        }

        if (instance.name.empty())
            instance.name = fileNameWithoutExtension(filePath);

        if (instance.dimension <= 0)
            throw std::runtime_error("TSPLIB file is missing a valid DIMENSION: " + filePath);

        if (maxDimensionForAdjMatrix > 0 && instance.dimension > maxDimensionForAdjMatrix)
        {
            throw std::runtime_error(
                "TSPLIB instance is too large for a dense adjacency matrix: "
                + instance.name
                + " has DIMENSION "
                + std::to_string(instance.dimension));
        }

        instance.edgeWeightType = normalizeTsplibToken(instance.edgeWeightType);
        instance.edgeWeightFormat = normalizeTsplibToken(instance.edgeWeightFormat);
        instance.optimalDist = knownOptimalDistance(instance.name);

        if (instance.edgeWeightType == "EXPLICIT")
        {
            instance.adjMat = buildExplicitMatrix(instance.dimension, instance.edgeWeightFormat, edgeWeights);
        }
        else
        {
            if ((int)points.size() != instance.dimension)
                throw std::runtime_error("TSPLIB coordinate count does not match DIMENSION in: " + filePath);

            instance.adjMat = buildCoordinateMatrix(instance.dimension, instance.edgeWeightType, points);
        }

        return instance;
    }

    static int knownOptimalDistance(const std::string& instanceName)
    {
        const auto& optima = knownOptima();
        const std::string key = toLower(fileNameWithoutExtension(instanceName));
        const auto it = optima.find(key);

        if (it == optima.end())
            return -1;

        return it->second;
    }

private:
    static void readHeaderLine(const std::string& line, TSPLibInstance& instance)
    {
        std::string key;
        std::string value;

        const size_t colon = line.find(':');
        if (colon != std::string::npos)
        {
            key = trim(line.substr(0, colon));
            value = trim(line.substr(colon + 1));
        }
        else
        {
            std::istringstream iss(line);
            iss >> key;
            std::getline(iss, value);
            value = trim(value);
        }

        key = normalizeHeaderKey(key);

        if (key == "NAME")
            instance.name = value;
        else if (key == "TYPE")
            instance.type = value;
        else if (key == "DIMENSION")
            instance.dimension = std::stoi(value);
        else if (key == "EDGEWEIGHTTYPE")
            instance.edgeWeightType = value;
        else if (key == "EDGEWEIGHTFORMAT")
            instance.edgeWeightFormat = value;
    }

    static void readPointLine(const std::string& line, int dimension, std::vector<Point>& points)
    {
        std::istringstream iss(line);
        int id = 0;
        double x = 0.0;
        double y = 0.0;

        if (!(iss >> id >> x >> y))
            return;

        if (id <= 0 || id > dimension)
            throw std::runtime_error("TSPLIB node id is outside DIMENSION.");

        if (points.empty())
            points.assign(dimension, Point());

        points[id - 1] = { x, y };
    }

    static void readIntLine(const std::string& line, std::vector<int>& values)
    {
        std::istringstream iss(line);
        int value = 0;

        while (iss >> value)
            values.push_back(value);
    }

    static std::vector<std::vector<int>> buildCoordinateMatrix(int n, const std::string& edgeWeightType, const std::vector<Point>& points)
    {
        std::vector<std::vector<int>> adj(n, std::vector<int>(n, 0));

        for (int i = 0; i < n; ++i)
        {
            for (int j = i + 1; j < n; ++j)
            {
                const int dist = coordinateDistance(edgeWeightType, points[i], points[j]);
                adj[i][j] = dist;
                adj[j][i] = dist;
            }
        }

        return adj;
    }

    static int coordinateDistance(const std::string& edgeWeightType, const Point& a, const Point& b)
    {
        const double dx = a.x - b.x;
        const double dy = a.y - b.y;

        if (edgeWeightType == "EUC_2D")
            return nint(std::sqrt(dx * dx + dy * dy));

        if (edgeWeightType == "CEIL_2D")
            return (int)std::ceil(std::sqrt(dx * dx + dy * dy));

        if (edgeWeightType == "ATT")
        {
            const double rij = std::sqrt((dx * dx + dy * dy) / 10.0);
            const int tij = nint(rij);
            return tij < rij ? tij + 1 : tij;
        }

        if (edgeWeightType == "GEO")
            return geoDistance(a, b);

        throw std::runtime_error("Unsupported TSPLIB EDGE_WEIGHT_TYPE: " + edgeWeightType);
    }

    static std::vector<std::vector<int>> buildExplicitMatrix(int n, const std::string& edgeWeightFormat, const std::vector<int>& values)
    {
        std::vector<std::vector<int>> adj(n, std::vector<int>(n, 0));
        size_t idx = 0;

        auto nextValue = [&]() -> int
        {
            if (idx >= values.size())
                throw std::runtime_error("EDGE_WEIGHT_SECTION ended before matrix was complete.");

            return values[idx++];
        };

        if (edgeWeightFormat.empty() || edgeWeightFormat == "FULL_MATRIX")
        {
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    adj[i][j] = nextValue();
        }
        else if (edgeWeightFormat == "UPPER_ROW")
        {
            for (int i = 0; i < n; ++i)
                for (int j = i + 1; j < n; ++j)
                    setSymmetric(adj, i, j, nextValue());
        }
        else if (edgeWeightFormat == "LOWER_ROW")
        {
            for (int i = 1; i < n; ++i)
                for (int j = 0; j < i; ++j)
                    setSymmetric(adj, i, j, nextValue());
        }
        else if (edgeWeightFormat == "UPPER_DIAG_ROW")
        {
            for (int i = 0; i < n; ++i)
                for (int j = i; j < n; ++j)
                    setSymmetric(adj, i, j, nextValue());
        }
        else if (edgeWeightFormat == "LOWER_DIAG_ROW")
        {
            for (int i = 0; i < n; ++i)
                for (int j = 0; j <= i; ++j)
                    setSymmetric(adj, i, j, nextValue());
        }
        else if (edgeWeightFormat == "UPPER_COL")
        {
            for (int j = 1; j < n; ++j)
                for (int i = 0; i < j; ++i)
                    setSymmetric(adj, i, j, nextValue());
        }
        else if (edgeWeightFormat == "LOWER_COL")
        {
            for (int j = 0; j < n - 1; ++j)
                for (int i = j + 1; i < n; ++i)
                    setSymmetric(adj, i, j, nextValue());
        }
        else if (edgeWeightFormat == "UPPER_DIAG_COL")
        {
            for (int j = 0; j < n; ++j)
                for (int i = 0; i <= j; ++i)
                    setSymmetric(adj, i, j, nextValue());
        }
        else if (edgeWeightFormat == "LOWER_DIAG_COL")
        {
            for (int j = 0; j < n; ++j)
                for (int i = j; i < n; ++i)
                    setSymmetric(adj, i, j, nextValue());
        }
        else
        {
            throw std::runtime_error("Unsupported TSPLIB EDGE_WEIGHT_FORMAT: " + edgeWeightFormat);
        }

        return adj;
    }

    static void setSymmetric(std::vector<std::vector<int>>& adj, int i, int j, int value)
    {
        adj[i][j] = value;
        adj[j][i] = value;
    }

    static int nint(double value)
    {
        return (int)(value + 0.5);
    }

    static int geoDistance(const Point& a, const Point& b)
    {
        constexpr double pi = 3.141592;
        constexpr double radius = 6378.388;

        const double latA = geoRadians(a.x, pi);
        const double lonA = geoRadians(a.y, pi);
        const double latB = geoRadians(b.x, pi);
        const double lonB = geoRadians(b.y, pi);

        const double q1 = std::cos(lonA - lonB);
        const double q2 = std::cos(latA - latB);
        const double q3 = std::cos(latA + latB);

        return (int)(radius * std::acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
    }

    static double geoRadians(double coordinate, double pi)
    {
        const int degrees = (int)coordinate;
        const double minutes = coordinate - degrees;
        return pi * (degrees + 5.0 * minutes / 3.0) / 180.0;
    }

    static std::string trim(const std::string& value)
    {
        size_t begin = 0;
        while (begin < value.size() && std::isspace((unsigned char)value[begin]))
            ++begin;

        size_t end = value.size();
        while (end > begin && std::isspace((unsigned char)value[end - 1]))
            --end;

        return value.substr(begin, end - begin);
    }

    static std::string toUpper(std::string value)
    {
        std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c)
        {
            return (char)std::toupper(c);
        });

        return value;
    }

    static std::string toLower(std::string value)
    {
        std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c)
        {
            return (char)std::tolower(c);
        });

        return value;
    }

    static std::string normalizeHeaderKey(std::string value)
    {
        value = toUpper(trim(value));
        value.erase(std::remove_if(value.begin(), value.end(), [](unsigned char c)
        {
            return std::isspace(c) || c == '_';
        }), value.end());

        return value;
    }

    static std::string normalizeTsplibToken(std::string value)
    {
        value = toUpper(trim(value));
        std::replace(value.begin(), value.end(), ' ', '_');
        return value;
    }

    static bool endsWith(const std::string& value, const std::string& suffix)
    {
        return value.size() >= suffix.size()
            && value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
    }

    static std::string fileNameWithoutExtension(const std::string& path)
    {
        const size_t slash = path.find_last_of("/\\");
        const size_t begin = slash == std::string::npos ? 0 : slash + 1;
        const size_t dot = path.find_last_of('.');

        if (dot == std::string::npos || dot < begin)
            return path.substr(begin);

        return path.substr(begin, dot - begin);
    }

    static const std::unordered_map<std::string, int>& knownOptima()
    {
        static const std::unordered_map<std::string, int> optima = {
            {"a280", 2579},
            {"ali535", 202310},
            {"att48", 10628},
            {"att532", 27686},
            {"bayg29", 1610},
            {"bays29", 2020},
            {"berlin52", 7542},
            {"bier127", 118282},
            {"brazil58", 25395},
            {"brg180", 1950},
            {"burma14", 3323},
            {"ch130", 6110},
            {"ch150", 6528},
            {"d198", 15780},
            {"d493", 35002},
            {"d657", 48912},
            {"d1291", 50801},
            {"d1655", 62128},
            {"dantzig42", 699},
            {"dsj1000", 18659688},
            {"eil51", 426},
            {"eil76", 538},
            {"fl417", 11861},
            {"fl1400", 20127},
            {"fnl4461", 182566},
            {"fri26", 937},
            {"gil262", 2378},
            {"gr17", 2085},
            {"gr21", 2707},
            {"gr24", 1272},
            {"gr48", 5046},
            {"gr96", 55209},
            {"gr120", 6942},
            {"gr137", 69853},
            {"gr202", 40160},
            {"gr229", 134602},
            {"gr431", 171414},
            {"gr666", 294358},
            {"hk48", 11461},
            {"kroa100", 21282},
            {"krob100", 22141},
            {"kroc100", 20749},
            {"krod100", 21294},
            {"kroe100", 22068},
            {"kroa150", 26524},
            {"krob150", 26130},
            {"kroa200", 29368},
            {"krob200", 29437},
            {"lin105", 14379},
            {"lin318", 42029},
            {"linhp318", 41345},
            {"nrw1379", 56638},
            {"p654", 34643},
            {"pa561", 2763},
            {"pcb442", 50778},
            {"pcb1173", 56892},
            {"pcb3038", 137694},
            {"pla7397", 23260728},
            {"pr76", 108159},
            {"pr107", 44303},
            {"pr124", 59030},
            {"pr136", 96772},
            {"pr144", 58537},
            {"pr152", 73682},
            {"pr226", 80369},
            {"pr264", 49135},
            {"pr299", 48191},
            {"pr439", 107217},
            {"pr1002", 259045},
            {"pr2392", 378032},
            {"rat99", 1211},
            {"rat195", 2323},
            {"rat575", 6773},
            {"rat783", 8806},
            {"rd100", 7910},
            {"rd400", 15281},
            {"rl1304", 252948},
            {"rl1323", 270199},
            {"rl1889", 316536},
            {"si175", 21407},
            {"si535", 48450},
            {"si1032", 92650},
            {"st70", 675},
            {"swiss42", 1273},
            {"ts225", 126643},
            {"tsp225", 3919},
            {"u159", 42080},
            {"u574", 36905},
            {"u724", 41910},
            {"u1060", 224094},
            {"u1432", 152970},
            {"u1817", 57201},
            {"u2152", 64253},
            {"u2319", 234256},
            {"ulysses16", 6859},
            {"ulysses22", 7013},
            {"vm1084", 239297},
            {"vm1748", 336556}
        };

        return optima;
    }
};
