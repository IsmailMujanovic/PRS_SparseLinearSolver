#include "../include/numeric.hpp"
#include <cmath>
#include <stdexcept>

NumericFactor::NumericFactor(std::size_t n) : diag(n, 0.0), lower(n) {}


NumericFactor numericFactorization(const Matrix& PAPt, const SymbolicFactorization& symbolic) {
    auto n = PAPt.n;

    NumericFactor L(n);
    L.lower.col_ptr = symbolic.col_ptr;
    L.lower.row_indices = symbolic.row_indices;
    L.lower.values.assign(symbolic.row_indices.size(), Complex(0.0, 0.0));

    struct Entry {
        std::size_t col;
        std::size_t pos;
    };
    std::vector<std::vector<Entry>> rowlist(n);
    for (std::size_t k = 0; k < n; k++) {
        for (auto idx = symbolic.col_ptr[k]; idx < symbolic.col_ptr[k + 1]; idx++) {
            rowlist[symbolic.row_indices[idx]].push_back({k, idx});
        }
    }

    std::vector<Complex> w(n, Complex(0.0, 0.0));
    for (std::size_t j = 0; j < n; j++) {
        for (auto idx = PAPt.col_ptr[j]; idx < PAPt.col_ptr[j + 1]; idx++) {
            w[PAPt.row_indices[idx]] = PAPt.values[idx];
        }

        for (const auto& entry : rowlist[j]) {
            auto conjLjk = std::conj(L.lower.values[entry.pos]);
            w[j] -= std::norm(L.lower.values[entry.pos]);

            for (auto idx = entry.pos + 1; idx < symbolic.col_ptr[entry.col + 1]; idx++) {
                w[symbolic.row_indices[idx]] -= L.lower.values[idx] * conjLjk;
            }
        }

        auto pivot = w[j].real();
        if (!std::isfinite(pivot) || pivot <= 0.0) {
            throw std::runtime_error("Matrix is not positive-definite.");
        }

        L.diag[j] = std::sqrt(pivot);
        w[j] = Complex(0.0, 0.0);
        for (auto idx = symbolic.col_ptr[j]; idx < symbolic.col_ptr[j + 1]; idx++) {
            auto r = symbolic.row_indices[idx];
            L.lower.values[idx] = w[r] / L.diag[j];
            w[r] = Complex(0.0, 0.0);
        }
    }

    return L;
}


std::vector<Complex> solve(const NumericFactor& L, const std::vector<std::size_t>& P, const std::vector<Complex>& b) {
    auto n = L.diag.size();
    const auto& lower = L.lower;

    std::vector<Complex> y(n);
    for (std::size_t i = 0; i < n; i++) {
        y[i] = b[P[i]];
    }

    for (std::size_t j = 0; j < n; j++) {
        y[j] /= L.diag[j];
        for (auto idx = lower.col_ptr[j]; idx < lower.col_ptr[j + 1]; idx++) {
            y[lower.row_indices[idx]] -= lower.values[idx] * y[j];
        }
    }

    for (auto j = n; j-- > 0;) {
        for (auto idx = lower.col_ptr[j]; idx < lower.col_ptr[j + 1]; idx++) {
            y[j] -= std::conj(lower.values[idx]) * y[lower.row_indices[idx]];
        }
        y[j] /= L.diag[j];
    }

    std::vector<Complex> x(n);
    for (std::size_t i = 0; i < n; i++) {
        x[P[i]] = y[i];
    }

    return x;
}
