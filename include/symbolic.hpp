#pragma once

#include <vector>
#include "matrix.hpp"

struct SymbolicFactorization {
    std::size_t n;
    std::vector<std::size_t> col_ptr;
    std::vector<std::size_t> row_indices;
    std::vector<std::size_t> parent;
};

SymbolicFactorization symbolicFactorization(const Matrix& PAPt);
