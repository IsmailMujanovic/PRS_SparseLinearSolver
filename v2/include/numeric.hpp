#pragma once

#include <vector>
#include "matrix.hpp"
#include "symbolic.hpp"

struct NumericFactor {
    std::vector<double> diag;
    Matrix lower;

    NumericFactor(std::size_t n);
};

NumericFactor numericFactorization(const Matrix& PAPt, const SymbolicFactorization& symbolic);

std::vector<Complex> solve(const NumericFactor& L, const std::vector<std::size_t>& P, const std::vector<Complex>& b);
