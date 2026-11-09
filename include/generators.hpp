#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <random>
#include <utility>
#include <vector>

#include "matrix.hpp"

inline Matrix randomHPD(std::size_t n, double density, unsigned seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> value(-1.0, 1.0);
    std::uniform_real_distribution<double> chance(0.0, 1.0);

    using DenseMatrix = std::vector<std::vector<Complex>>;
    DenseMatrix H(n, std::vector<Complex>(n));

    for (std::size_t i = 0; i < n; i++) {
        for (std::size_t j = 0; j < i; j++) {
            if (chance(gen) < density) {
                H[i][j] = Complex(value(gen), value(gen));
                H[j][i] = std::conj(H[i][j]);
            }
        }
    }

    for (std::size_t i = 0; i < n; i++) {
        double sum = 1.0;
        for (std::size_t j = 0; j < n; j++) {
            sum += std::abs(H[i][j]);
        }
        H[i][i] = sum;
    }

    Matrix A(n);
    for (std::size_t c = 0; c < n; c++) {
        for (std::size_t r = 0; r < n; r++) {
            if (H[r][c] != Complex(0.0, 0.0)) {
                A.row_indices.push_back(r);
                A.values.push_back(H[r][c]);
            }
        }
        A.col_ptr[c + 1] = A.row_indices.size();
    }

    return A;
}

inline Matrix tridiagonalHPD(std::size_t n) {
    Matrix A(n);
    std::vector<std::vector<std::pair<std::size_t, Complex>>> columns(n);

    for (std::size_t i = 0; i < n; i++) {
        columns[i].push_back({i, Complex(4.0, 0.0)});
        if (i > 0) {
            columns[i].push_back({i - 1, Complex(1.0, 1.0)});
        }
        if (i + 1 < n) {
            columns[i].push_back({i + 1, Complex(1.0, -1.0)});
        }
    }

    for (std::size_t c = 0; c < n; c++) {
        std::sort(columns[c].begin(), columns[c].end(), [](const auto& a, const auto& b) {
            return a.first < b.first;
        });
        for (const auto& entry : columns[c]) {
            A.row_indices.push_back(entry.first);
            A.values.push_back(entry.second);
        }
        A.col_ptr[c + 1] = A.row_indices.size();
    }

    return A;
}

inline Matrix complexLaplacian2D(std::size_t k) {
    std::size_t n = k * k;
    std::vector<std::vector<std::pair<std::size_t, Complex>>> columns(n);

    auto index = [k](std::size_t i, std::size_t j) {
        return i * k + j;
    };

    auto addEdge = [&columns](std::size_t a, std::size_t b, double theta) {
        Complex weight = -std::polar(1.0, theta);
        if (a < b) {
            columns[b].push_back({a, weight});
            columns[a].push_back({b, std::conj(weight)});
        } else {
            columns[a].push_back({b, weight});
            columns[b].push_back({a, std::conj(weight)});
        }
    };

    for (std::size_t i = 0; i < k; i++) {
        for (std::size_t j = 0; j < k; j++) {
            std::size_t node = index(i, j);
            if (j + 1 < k) {
                double theta = 0.37 + 0.05 * static_cast<double>((i + j) % 5);
                addEdge(node, index(i, j + 1), theta);
            }
            if (i + 1 < k) {
                double theta = -0.81 + 0.07 * static_cast<double>((i + 2 * j) % 5);
                addEdge(node, index(i + 1, j), theta);
            }
        }
    }

    for (std::size_t c = 0; c < n; c++) {
        double diagonal = 1.0;
        for (const auto& entry : columns[c]) {
            diagonal += std::abs(entry.second);
        }
        columns[c].push_back({c, Complex(diagonal, 0.0)});
        std::sort(columns[c].begin(), columns[c].end(), [](const auto& a, const auto& b) {
            return a.first < b.first;
        });
    }

    Matrix A(n);
    for (std::size_t c = 0; c < n; c++) {
        for (const auto& entry : columns[c]) {
            A.row_indices.push_back(entry.first);
            A.values.push_back(entry.second);
        }
        A.col_ptr[c + 1] = A.row_indices.size();
    }

    return A;
}
