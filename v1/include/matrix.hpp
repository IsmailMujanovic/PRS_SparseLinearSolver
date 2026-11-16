#pragma once

#include <vector>
#include <complex>

using Complex = std::complex<double>;

struct Matrix {
    std::size_t n;
    std::vector<Complex> values;
    std::vector<std::size_t> row_indices;
    std::vector<std::size_t> col_ptr;
    
    Matrix(std::size_t n);
};