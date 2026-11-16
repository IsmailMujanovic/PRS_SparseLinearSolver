#include "../include/matrix.hpp"
#include "../include/amd.hpp"
#include "../include/permute.hpp"
#include "../include/symbolic.hpp"
#include "../include/numeric.hpp"
#include "../include/validate.hpp"
#include "../include/generators.hpp"
#include <algorithm>
#include <chrono>
#include <iostream>

using Dense = std::vector<std::vector<Complex>>;

#ifdef PRS_ENABLE_PROFILE
struct PhaseTimings {
    double validation = 0.0;
    double ordering = 0.0;
    double permutation = 0.0;
    double symbolic = 0.0;
    double numeric = 0.0;
    double solve = 0.0;
};

double elapsedSeconds(const std::chrono::steady_clock::time_point& start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}
#endif

Matrix toCSC(const Dense& H) {
    auto n = H.size();

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
#ifdef PRS_ENABLE_PROFILE
    PhaseTimings timings;
    auto phaseStart = std::chrono::steady_clock::now();
#endif
    auto validation = validateHPD(A);
#ifdef PRS_ENABLE_PROFILE
    timings.validation = elapsedSeconds(phaseStart);
#endif
    if (!validation.ok()) {
        std::cout << name << ": odbijena matrica (" << validation.message << ")\n";
        return false;
    }

    auto n = A.n;

    // Build RHS vector 'b' BEFORE starting the profiler so it's always defined
    // and doesn't pollute solver profile samples
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

#ifdef PRS_ENABLE_PROFILE
    phaseStart = std::chrono::steady_clock::now();
#endif
    auto P = AMD(A);
#ifdef PRS_ENABLE_PROFILE
    timings.ordering = elapsedSeconds(phaseStart);
    phaseStart = std::chrono::steady_clock::now();
#endif
    auto PAPt = permute(A, P);
#ifdef PRS_ENABLE_PROFILE
    timings.permutation = elapsedSeconds(phaseStart);
    phaseStart = std::chrono::steady_clock::now();
#endif
    auto symbolic = symbolicFactorization(PAPt);
#ifdef PRS_ENABLE_PROFILE
    timings.symbolic = elapsedSeconds(phaseStart);
    phaseStart = std::chrono::steady_clock::now();
#endif
    auto L = numericFactorization(PAPt, symbolic);
#ifdef PRS_ENABLE_PROFILE
    timings.numeric = elapsedSeconds(phaseStart);
    phaseStart = std::chrono::steady_clock::now();
#endif
    auto x = solve(L, P, b);
#ifdef PRS_ENABLE_PROFILE
    timings.solve = elapsedSeconds(phaseStart);
#endif

#ifndef PRS_ENABLE_PROFILE
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

    double solveError = 0.0;
    for (std::size_t i = 0; i < n; i++) {
        solveError = std::max(solveError, std::abs(x[i] - xKnown[i]));
    }

    std::cout << name << ": n=" << n << ", nnz(A)=" << A.values.size() << ", nnz(L)=" << L.lower.values.size() + n
              << ", |PAP^H-LL^H|=" << factorError << ", |x-xKnown|=" << solveError << "\n";

    return factorError < 1e-10 && solveError < 1e-10;
#else
    // Prevent compiler from optimizing away 'x' when error check is disabled
    std::cout << name << ": n=" << n << ", nnz(A)=" << A.values.size() << ", nnz(L)=" << L.lower.values.size() + n
              << ", x[0]=" << x[0] << "\n";
    std::cout << "  profile [ms]: validation=" << timings.validation * 1000.0
              << ", ordering=" << timings.ordering * 1000.0
              << ", permutation=" << timings.permutation * 1000.0
              << ", symbolic=" << timings.symbolic * 1000.0
              << ", numeric=" << timings.numeric * 1000.0
              << ", solve=" << timings.solve * 1000.0 << "\n";
    return true;
#endif
}

bool expectRejected(const Matrix& A, const std::string& name) {
    auto validation = validateHPD(A);
    if (validation.ok()) {
        std::cout << name << ": GREŠKA - matrica je prihvaćena, a trebala je biti odbijena\n";
        return false;
    }

    std::cout << name << ": odbijena kako je očekivano (" << validation.message << ")\n";
    return true;
}

int main() {
    auto ok = check(randomHPD(5000, 0.02, 2), "random 5000");

#ifndef PRS_ENABLE_PROFILE
    Matrix small(4);
    small.values = {
        {4.0, 0.0}, {1.0, -1.0},
        {1.0, 1.0}, {5.0, 0.0}, {0.0, -2.0},
                    {0.0, 2.0}, {4.0, 0.0}, {1.0, 0.0},
                                {1.0, 0.0}, {3.0, 0.0}};
    small.row_indices = {0, 1, 0, 1, 2, 1, 2, 3, 2, 3};
    small.col_ptr = {0, 2, 5, 8, 10};

    Dense nonHermitian = {
        {{1.0, 0.0}, {2.0, 1.0}},
        {{3.0, -1.0}, {4.0, 0.0}}};
    Dense indefinite = {
        {1.0, 2.0},
        {2.0, 1.0}};
    Dense nonRealDiagonal = {
        {{1.0, 1.0}, {0.0, 0.0}},
        {0.0, 2.0}};
    Dense complexHpd = {
        {2.0, {1.0, 1.0}},
        {{1.0, -1.0}, 3.0}};

    ok = check(small, "small 4x4") && ok;
    ok = check(randomHPD(60, 0.08, 1), "random 60") && ok;
    ok = expectRejected(toCSC(nonHermitian), "non-Hermitian") && ok;
    ok = expectRejected(toCSC(indefinite), "indefinite") && ok;
    ok = expectRejected(toCSC(nonRealDiagonal), "non-real diagonal") && ok;
    ok = check(toCSC(complexHpd), "complex HPD 2x2") && ok;
#endif

    std::cout << (ok ? "ALL OK" : "FAILED") << "\n";
    return ok ? 0 : 1;
}
