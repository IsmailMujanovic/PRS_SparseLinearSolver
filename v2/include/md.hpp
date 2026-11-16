#pragma once

#include <vector>
#include <set>
#include "matrix.hpp"

class MD_Graph {
private:
    std::size_t n;
    std::vector<std::size_t> degrees;
    std::vector<std::set<std::size_t>> edges;
    std::vector<bool> active;

public:
    MD_Graph(const Matrix& matrix);
    std::size_t findMinimumDegree() const;
    void eliminateVertex(std::size_t vertex);
    bool isActive(std::size_t vertex) const;
};

std::vector<std::size_t> MD(const Matrix& matrix);