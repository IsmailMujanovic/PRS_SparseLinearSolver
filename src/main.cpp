#include "../include/matrix.hpp"
#include "../include/md.hpp"
#include "../include/permute.hpp"
#include "../include/symbolic.hpp"
#include "../include/numeric.hpp"
#include <algorithm>
#include <iostream>
#include <random>

using Dense = std::vector<std::vector<Complex>>;

Matrix randomHPD(std::size_t n, double density, unsigned seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> value(-1.0, 1.0);
    std::uniform_real_distribution<double> chance(0.0, 1.0);

    Dense H(n, std::vector<Complex>(n));
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

bool check(const Matrix& A, const std::string& name) {
    auto n = A.n;
    auto P = MD(A);
    auto PAPt = permute(A, P);
    auto symbolic = symbolicFactorization(PAPt);
    auto L = numericFactorization(PAPt, symbolic);

    Dense dense(n, std::vector<Complex>(n));
    for (std::size_t c = 0; c < n; c++) {
        for (auto idx = PAPt.col_ptr[c]; idx < PAPt.col_ptr[c + 1]; idx++) {
            auto r = PAPt.row_indices[idx];
            dense[r][c] = PAPt.values[idx];
            dense[c][r] = std::conj(PAPt.values[idx]);
        }
    }

    Dense denseL(n, std::vector<Complex>(n));
    for (std::size_t c = 0; c < n; c++) {
        denseL[c][c] = L.diag[c];
        for (auto idx = L.lower.col_ptr[c]; idx < L.lower.col_ptr[c + 1]; idx++) {
            denseL[L.lower.row_indices[idx]][c] = L.lower.values[idx];
        }
    }

    double factorError = 0.0;
    for (std::size_t i = 0; i < n; i++) {
        for (std::size_t j = 0; j < n; j++) {
            Complex sum(0.0, 0.0);
            for (std::size_t k = 0; k < n; k++) {
                sum += denseL[i][k] * std::conj(denseL[j][k]);
            }
            factorError = std::max(factorError, std::abs(sum - dense[i][j]));
        }
    }

    std::vector<Complex> xKnown(n);
    for (std::size_t i = 0; i < n; i++) {
        xKnown[i] = Complex(i + 1.0, 0.5 * i - 1.0);
    }

    std::vector<Complex> b(n);
    for (std::size_t c = 0; c < n; c++) {
        for (auto idx = A.col_ptr[c]; idx < A.col_ptr[c + 1]; idx++) {
            b[A.row_indices[idx]] += A.values[idx] * xKnown[c];
        }
    }

    auto x = solve(L, P, b);
    double solveError = 0.0;
    for (std::size_t i = 0; i < n; i++) {
        solveError = std::max(solveError, std::abs(x[i] - xKnown[i]));
    }

    std::cout << name << ": n=" << n << ", nnz(A)=" << A.values.size() << ", nnz(L)=" << L.lower.values.size() + n
              << ", |PAP^H-LL^H|=" << factorError << ", |x-xKnown|=" << solveError << "\n";

    return factorError < 1e-10 && solveError < 1e-10;
}

int main() {
    Matrix small(4);
    small.values = {
        {4.0, 0.0}, {1.0, -1.0},
        {1.0, 1.0}, {5.0, 0.0}, {0.0, -2.0},
                    {0.0, 2.0}, {4.0, 0.0}, {1.0, 0.0},
                                {1.0, 0.0}, {3.0, 0.0}};
    small.row_indices = {0, 1, 0, 1, 2, 1, 2, 3, 2, 3};
    small.col_ptr = {0, 2, 5, 8, 10};

    auto ok = check(small, "tridiagonal");
    ok = check(randomHPD(60, 0.08, 1), "random 60") && ok;
    ok = check(randomHPD(300, 0.02, 2), "random 300") && ok;

    std::cout << (ok ? "ALL OK" : "FAILED") << "\n";
    return ok ? 0 : 1;
}
