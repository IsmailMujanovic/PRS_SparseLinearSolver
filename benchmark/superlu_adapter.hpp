#pragma once

#include <cstddef>
#include <vector>

#include "../include/matrix.hpp"

enum class SuperLUOrderingMethod {
    MmdAtPlusA,
    Colamd
};

enum class SuperLUMode {
    Default,
    Hpd,
    SameOrder
};

struct SuperLUOrderingResult {
    std::vector<std::size_t> permutation;
    double seconds = 0.0;
};

struct SuperLUOptions {
    SuperLUMode mode = SuperLUMode::Hpd;
    const std::vector<std::size_t>* ordering = nullptr;
};

struct SuperLUPhaseTimes {
    double ordering = 0.0;
    double symbolic = 0.0;
    double numeric = 0.0;
    double solve = 0.0;
    double total = 0.0;
    double etree = 0.0;
    double fact_first = 0.0;
    double fact_reuse = 0.0;
    double wall = 0.0;
};

struct SuperLUResult {
    std::vector<Complex> solution;
    SuperLUPhaseTimes times;
    std::size_t nnz_l = 0;
    std::size_t nnz_u = 0;
    std::size_t nnz_factor = 0;
    double factor_bytes = 0.0;
    double total_bytes = 0.0;
    double flops = 0.0;
    int info = 0;
    bool split_valid = false;
};

SuperLUOrderingResult superluOrdering(const Matrix& A, SuperLUOrderingMethod method);

SuperLUResult runSuperLU(const Matrix& A, const std::vector<Complex>& b, const SuperLUOptions& options);
