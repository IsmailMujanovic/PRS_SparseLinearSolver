#include "../include/matrix.hpp"

Matrix::Matrix(std::size_t n) : n(n), values(), row_indices(), col_ptr(n + 1) {}