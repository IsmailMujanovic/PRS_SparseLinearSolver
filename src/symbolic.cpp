#include "../include/symbolic.hpp"
#include <algorithm>
#include <limits>

SymbolicFactorization symbolicFactorization(const Matrix& PAPt) {
    auto n = PAPt.n;

    SymbolicFactorization symbolic;
    symbolic.n = n;
    symbolic.col_ptr.assign(n + 1, 0);
    symbolic.parent.assign(n, n);

    std::vector<std::size_t> mark(n, std::numeric_limits<std::size_t>::max());
    std::vector<std::vector<std::size_t>> children(n);

    for (std::size_t j = 0; j < n; j++) {
        auto start = symbolic.row_indices.size();
        auto minRow = n;

        auto add = [&](std::size_t r) {
            if (r > j && mark[r] != j) {
                mark[r] = j;
                symbolic.row_indices.push_back(r);
                minRow = std::min(minRow, r);
            }
        };

        for (auto idx = PAPt.col_ptr[j]; idx < PAPt.col_ptr[j + 1]; idx++) {
            add(PAPt.row_indices[idx]);
        }

        for (auto child : children[j]) {
            for (auto idx = symbolic.col_ptr[child]; idx < symbolic.col_ptr[child + 1]; idx++) {
                add(symbolic.row_indices[idx]);
            }
        }

        std::sort(symbolic.row_indices.begin() + start, symbolic.row_indices.end());
        symbolic.col_ptr[j + 1] = symbolic.row_indices.size();

        if (minRow < n) {
            symbolic.parent[j] = minRow;
            children[minRow].push_back(j);
        }
    }

    return symbolic;
}
