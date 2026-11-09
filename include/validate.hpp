#pragma once

#include <string>
#include "matrix.hpp"

struct ValidationResult {
    bool hermitian = false;
    bool positiveDefinite = false;
    std::string message;

    bool ok() const;
};

ValidationResult validateHPD(const Matrix& A, double tolerance = 1e-12);
