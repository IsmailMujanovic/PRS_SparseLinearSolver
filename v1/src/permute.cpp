#include "../include/permute.hpp"
#include <algorithm>
#include <utility>

Matrix permute(const Matrix& A, const std::vector<std::size_t>& P) {
    auto n = A.n;

    std::vector<std::size_t> invP(n);
    for (std::size_t i = 0; i < n; i++) {
        invP[P[i]] = i;
    }

    std::vector<std::vector<std::pair<std::size_t, Complex>>> columns(n);
    for (std::size_t c = 0; c < n; c++) {
        for (auto idx = A.col_ptr[c]; idx < A.col_ptr[c + 1]; idx++) {
            auto r = A.row_indices[idx];
            if (r < c) {
                continue;
            }

            auto newR = invP[r];
            auto newC = invP[c];
            if (newR >= newC) {
                columns[newC].emplace_back(newR, A.values[idx]);
            } else {
                columns[newR].emplace_back(newC, std::conj(A.values[idx]));
            }
        }
    }

    Matrix result(n);
    for (std::size_t j = 0; j < n; j++) {
        std::sort(columns[j].begin(), columns[j].end(), [](const auto& a, const auto& b) {
            return a.first < b.first;
        });

        for (const auto& [row, value] : columns[j]) {
            result.row_indices.push_back(row);
            result.values.push_back(value);
        }
        result.col_ptr[j + 1] = result.row_indices.size();
    }

    return result;
}
