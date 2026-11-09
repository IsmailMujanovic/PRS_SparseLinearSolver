#include "../include/validate.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>
#include <utility>
#include <vector>

namespace {

std::string formatComplex(const Complex& z) {
    std::ostringstream out;
    out << "(" << z.real() << ", " << z.imag() << ")";
    return out.str();
}

std::string formatEntry(std::size_t r, std::size_t c) {
    return "a[" + std::to_string(r) + "][" + std::to_string(c) + "]";
}

bool sameValue(const Complex& a, const Complex& b, double tolerance) {
    double scale = std::abs(a) + std::abs(b) + 1.0;
    return std::abs(a - b) <= tolerance * scale;
}

}

bool ValidationResult::ok() const {
    return hermitian && positiveDefinite;
}


ValidationResult validateHPD(const Matrix& A, double tolerance) {
    ValidationResult result;

    if (A.col_ptr.size() != A.n + 1) {
        result.message = "neispravna CSC struktura: col_ptr mora imati n+1 elemenata";
        return result;
    }
    if (A.values.size() != A.row_indices.size()) {
        result.message = "neispravna CSC struktura: broj vrijednosti i indeksa redova se razlikuje";
        return result;
    }
    if (A.col_ptr.front() != 0 || A.col_ptr.back() != A.values.size()) {
        result.message = "neispravna CSC struktura: col_ptr ne pokriva sve vrijednosti";
        return result;
    }
    for (std::size_t c = 0; c < A.n; c++) {
        if (A.col_ptr[c] > A.col_ptr[c + 1]) {
            result.message = "neispravna CSC struktura: col_ptr nije monotono rastući";
            return result;
        }
    }

    std::map<std::pair<std::size_t, std::size_t>, Complex> entries;
    for (std::size_t c = 0; c < A.n; c++) {
        for (auto idx = A.col_ptr[c]; idx < A.col_ptr[c + 1]; idx++) {
            auto r = A.row_indices[idx];
            if (r >= A.n) {
                result.message = "neispravna CSC struktura: indeks reda je van opsega";
                return result;
            }
            entries[{r, c}] += A.values[idx];
        }
    }

    for (const auto& [position, value] : entries) {
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) {
            result.message = formatEntry(position.first, position.second) + " nije konačan";
            return result;
        }
    }

    for (const auto& [position, value] : entries) {
        const auto& [r, c] = position;
        if (r == c) {
            if (std::abs(value.imag()) > tolerance * (std::abs(value) + 1.0)) {
                result.message = "dijagonalni element " + formatEntry(r, c) + " = " +
                                 formatComplex(value) + " nije realan";
                return result;
            }
            continue;
        }

        auto mirrorIt = entries.find({c, r});
        Complex mirror = mirrorIt == entries.end() ? Complex(0.0, 0.0) : mirrorIt->second;
        if (!sameValue(value, std::conj(mirror), tolerance)) {
            result.message = "matrica nije Hermitijska: " + formatEntry(r, c) + " = " +
                             formatComplex(value) + ", a " + formatEntry(c, r) + " = " +
                             formatComplex(mirror);
            return result;
        }
    }
    result.hermitian = true;

    std::vector<std::vector<Complex>> dense(A.n, std::vector<Complex>(A.n, Complex(0.0, 0.0)));
    for (const auto& [position, value] : entries) {
        dense[position.first][position.second] = value;
    }

    double scale = 0.0;
    for (std::size_t i = 0; i < A.n; i++) {
        scale = std::max(scale, std::abs(dense[i][i]));
    }
    if (scale == 0.0) {
        scale = 1.0;
    }

    std::vector<std::vector<Complex>> L(A.n, std::vector<Complex>(A.n, Complex(0.0, 0.0)));
    for (std::size_t k = 0; k < A.n; k++) {
        Complex pivotSum = dense[k][k];
        for (std::size_t s = 0; s < k; s++) {
            pivotSum -= L[k][s] * std::conj(L[k][s]);
        }

        double pivot = pivotSum.real();
        if (!std::isfinite(pivot) || pivot <= tolerance * scale) {
            result.message = "matrica nije pozitivno definitna: pivot na indeksu " +
                             std::to_string(k) + " je " + std::to_string(pivot);
            return result;
        }

        L[k][k] = std::sqrt(pivot);
        for (std::size_t i = k + 1; i < A.n; i++) {
            Complex sum = dense[i][k];
            for (std::size_t s = 0; s < k; s++) {
                sum -= L[i][s] * std::conj(L[k][s]);
            }
            L[i][k] = sum / L[k][k];
        }
    }

    result.positiveDefinite = true;
    return result;
}
