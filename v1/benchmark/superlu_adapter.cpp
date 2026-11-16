#include "superlu_adapter.hpp"

#include <slu_zdefs.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

template <typename T>
struct SluDeleter {
    void operator()(T* pointer) const {
        if (pointer != nullptr) {
            superlu_free(pointer);
        }
    }
};

template <typename T>
using SluPtr = std::unique_ptr<T, SluDeleter<T>>;

struct MatrixGuard {
    SuperMatrix* matrix;
    bool compressedColumn;

    ~MatrixGuard() {
        if (matrix == nullptr) {
            return;
        }
        if (compressedColumn) {
            Destroy_CompCol_Matrix(matrix);
        } else {
            Destroy_SuperMatrix_Store(matrix);
        }
    }
};

struct FactorGuard {
    SuperMatrix* lower;
    SuperMatrix* upper;
    bool active;

    ~FactorGuard() {
        if (!active) {
            return;
        }
        Destroy_SuperNode_Matrix(lower);
        Destroy_CompCol_Matrix(upper);
    }
};

struct StatGuard {
    SuperLUStat_t* stat;

    ~StatGuard() {
        StatFree(stat);
    }
};

int checkedSize(std::size_t value, const char* what) {
    if (value > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error(std::string(what) + " does not fit into SuperLU int indices");
    }
    return static_cast<int>(value);
}

int_t checkedNonzeros(std::size_t value) {
    if (value > static_cast<std::size_t>(std::numeric_limits<int_t>::max())) {
        throw std::runtime_error("matrix nonzeros do not fit into SuperLU int_t");
    }
    return static_cast<int_t>(value);
}

SuperMatrix createCompressedColumnMatrix(const Matrix& A) {
    int n = checkedSize(A.n, "matrix dimension");
    int_t nnz = checkedNonzeros(A.values.size());

    std::size_t valueCount = nnz > 0 ? static_cast<std::size_t>(nnz) : 1;
    SluPtr<doublecomplex> values(doublecomplexMalloc(valueCount));
    SluPtr<int_t> rowIndices(intMalloc(nnz > 0 ? nnz : 1));
    SluPtr<int_t> columnPointers(intMalloc(static_cast<int_t>(n) + 1));
    if (!values || !rowIndices || !columnPointers) {
        throw std::runtime_error("SuperLU allocation failed while creating the matrix");
    }

    for (std::size_t index = 0; index < A.values.size(); index++) {
        values.get()[index].r = A.values[index].real();
        values.get()[index].i = A.values[index].imag();
        rowIndices.get()[index] = static_cast<int_t>(A.row_indices[index]);
    }
    for (int column = 0; column <= n; column++) {
        columnPointers.get()[column] = static_cast<int_t>(A.col_ptr[static_cast<std::size_t>(column)]);
    }

    SuperMatrix matrix;
    zCreate_CompCol_Matrix(&matrix, n, n, nnz, values.release(), rowIndices.release(),
                           columnPointers.release(), SLU_NC, SLU_Z, SLU_GE);
    return matrix;
}

void validatePermutation(const std::vector<std::size_t>& permutation, int n) {
    if (permutation.size() != static_cast<std::size_t>(n)) {
        throw std::runtime_error("ordering size does not match the matrix dimension");
    }
    std::vector<char> seen(static_cast<std::size_t>(n), 0);
    for (std::size_t position : permutation) {
        if (position >= static_cast<std::size_t>(n) || seen[position] != 0) {
            throw std::runtime_error("ordering is not a valid permutation");
        }
        seen[position] = 1;
    }
}

}

SuperLUOrderingResult superluOrdering(const Matrix& A, SuperLUOrderingMethod method) {
    int n = checkedSize(A.n, "matrix dimension");
    if (n == 0) {
        return SuperLUOrderingResult{};
    }

    SuperMatrix matrix = createCompressedColumnMatrix(A);
    MatrixGuard matrixGuard{&matrix, true};

    SluPtr<int> permutation(int32Malloc(n));
    if (!permutation) {
        throw std::runtime_error("SuperLU allocation failed for the ordering");
    }

    int specification = method == SuperLUOrderingMethod::MmdAtPlusA ? MMD_AT_PLUS_A : COLAMD;
    auto start = std::chrono::steady_clock::now();
    get_perm_c(specification, &matrix, permutation.get());
    auto stop = std::chrono::steady_clock::now();

    SuperLUOrderingResult result;
    result.seconds = std::chrono::duration<double>(stop - start).count();
    result.permutation.resize(static_cast<std::size_t>(n));
    std::vector<char> seen(static_cast<std::size_t>(n), 0);
    for (int column = 0; column < n; column++) {
        int position = permutation.get()[column];
        if (position < 0 || position >= n || seen[position] != 0) {
            throw std::runtime_error("SuperLU returned an invalid column permutation");
        }
        seen[position] = 1;
        result.permutation[static_cast<std::size_t>(position)] = static_cast<std::size_t>(column);
    }
    return result;
}

SuperLUResult runSuperLU(const Matrix& A, const std::vector<Complex>& b, const SuperLUOptions& options) {
    int n = checkedSize(A.n, "matrix dimension");
    if (b.size() != A.n) {
        throw std::runtime_error("right-hand side size does not match the matrix dimension");
    }
    if (n == 0) {
        return SuperLUResult{};
    }

    SuperMatrix matrix = createCompressedColumnMatrix(A);
    MatrixGuard matrixGuard{&matrix, true};

    SluPtr<int> permC(int32Malloc(n));
    SluPtr<int> permR(int32Malloc(n));
    SluPtr<int> etree(int32Malloc(n));
    SluPtr<double> rowScale(doubleMalloc(n));
    SluPtr<double> columnScale(doubleMalloc(n));
    SluPtr<double> forwardError(doubleMalloc(1));
    SluPtr<double> backwardError(doubleMalloc(1));
    SluPtr<doublecomplex> rhs(doublecomplexMalloc(static_cast<std::size_t>(n)));
    SluPtr<doublecomplex> solution(doublecomplexMalloc(static_cast<std::size_t>(n)));
    if (!permC || !permR || !etree || !rowScale || !columnScale || !forwardError ||
        !backwardError || !rhs || !solution) {
        throw std::runtime_error("SuperLU allocation failed for the solver workspace");
    }

    for (int index = 0; index < n; index++) {
        rhs.get()[index].r = b[static_cast<std::size_t>(index)].real();
        rhs.get()[index].i = b[static_cast<std::size_t>(index)].imag();
        solution.get()[index].r = 0.0;
        solution.get()[index].i = 0.0;
    }

    SuperMatrix rightHandSide;
    SuperMatrix solutionMatrix;
    zCreate_Dense_Matrix(&rightHandSide, n, 1, rhs.get(), n, SLU_DN, SLU_Z, SLU_GE);
    MatrixGuard rightHandSideGuard{&rightHandSide, false};
    zCreate_Dense_Matrix(&solutionMatrix, n, 1, solution.get(), n, SLU_DN, SLU_Z, SLU_GE);
    MatrixGuard solutionMatrixGuard{&solutionMatrix, false};

    superlu_options_t solverOptions;
    set_default_options(&solverOptions);
    solverOptions.Equil = NO;
    solverOptions.IterRefine = NOREFINE;
    solverOptions.PrintStat = NO;
    solverOptions.PivotGrowth = NO;
    solverOptions.ConditionNumber = NO;
    solverOptions.Fact = DOFACT;

    if (options.mode == SuperLUMode::Default) {
        solverOptions.ColPerm = COLAMD;
        solverOptions.DiagPivotThresh = 1.0;
        solverOptions.SymmetricMode = NO;
    } else {
        solverOptions.DiagPivotThresh = 0.0;
        solverOptions.SymmetricMode = YES;
        if (options.mode == SuperLUMode::Hpd) {
            solverOptions.ColPerm = MMD_AT_PLUS_A;
        } else {
            solverOptions.ColPerm = MY_PERMC;
            if (options.ordering == nullptr) {
                throw std::runtime_error("same-ordering mode requires a permutation");
            }
            validatePermutation(*options.ordering, n);
            for (std::size_t position = 0; position < options.ordering->size(); position++) {
                permC.get()[(*options.ordering)[position]] = static_cast<int>(position);
            }
        }
    }

    SuperMatrix lower;
    SuperMatrix upper;
    FactorGuard factorGuard{&lower, &upper, false};

    SuperLUStat_t firstStat;
    StatInit(&firstStat);
    StatGuard firstStatGuard{&firstStat};

    GlobalLU_t globalLU;
    mem_usage_t memoryUsage;
    double reciprocalPivotGrowth = 0.0;
    double reciprocalCondition = 0.0;
    char equilibration = 'N';
    int_t info = 0;

    auto start = std::chrono::steady_clock::now();
    zgssvx(&solverOptions, &matrix, permC.get(), permR.get(), etree.get(), &equilibration,
           rowScale.get(), columnScale.get(), &lower, &upper, nullptr, 0, &rightHandSide,
           &solutionMatrix, &reciprocalPivotGrowth, &reciprocalCondition, forwardError.get(),
           backwardError.get(), &globalLU, &memoryUsage, &firstStat, &info);
    auto stop = std::chrono::steady_clock::now();

    if (info >= 0 && info <= n + 1) {
        factorGuard.active = true;
    }
    if (info < 0) {
        throw std::runtime_error("SuperLU reported an invalid argument, info=" + std::to_string(info));
    }
    if (info > n + 1) {
        throw std::runtime_error("SuperLU memory allocation failed, info=" + std::to_string(info));
    }
    if (info > 0 && info <= n) {
        throw std::runtime_error("SuperLU detected a singular factor, info=" + std::to_string(info));
    }

    SuperLUResult result;
    result.info = static_cast<int>(info);
    result.times.ordering = firstStat.utime[COLPERM];
    result.times.etree = firstStat.utime[ETREE];
    result.times.fact_first = firstStat.utime[FACT];
    result.times.solve = firstStat.utime[SOLVE];
    result.times.wall = std::chrono::duration<double>(stop - start).count();
    result.flops = static_cast<double>(firstStat.ops[FACT]);
    result.factor_bytes = static_cast<double>(memoryUsage.for_lu);
    result.total_bytes = static_cast<double>(memoryUsage.total_needed);
    result.split_valid = options.mode != SuperLUMode::Default;

    auto* lowerStore = static_cast<SCformat*>(lower.Store);
    auto* upperStore = static_cast<NCformat*>(upper.Store);
    result.nnz_l = static_cast<std::size_t>(lowerStore->nnz);
    result.nnz_u = static_cast<std::size_t>(upperStore->nnz);
    result.nnz_factor = result.nnz_l + result.nnz_u - static_cast<std::size_t>(n);

    result.solution.resize(static_cast<std::size_t>(n));
    auto* solutionValues = static_cast<doublecomplex*>(static_cast<DNformat*>(solutionMatrix.Store)->nzval);
    for (int index = 0; index < n; index++) {
        result.solution[static_cast<std::size_t>(index)] =
            Complex(solutionValues[index].r, solutionValues[index].i);
    }

    SuperLUStat_t secondStat;
    StatInit(&secondStat);
    StatGuard secondStatGuard{&secondStat};

    solverOptions.Fact = SamePattern_SameRowPerm;
    rightHandSide.ncol = 0;
    solutionMatrix.ncol = 0;
    mem_usage_t reuseMemoryUsage;
    int_t reuseInfo = 0;

    zgssvx(&solverOptions, &matrix, permC.get(), permR.get(), etree.get(), &equilibration,
           rowScale.get(), columnScale.get(), &lower, &upper, nullptr, 0, &rightHandSide,
           &solutionMatrix, &reciprocalPivotGrowth, &reciprocalCondition, forwardError.get(),
           backwardError.get(), &globalLU, &reuseMemoryUsage, &secondStat, &reuseInfo);

    if (reuseInfo < 0) {
        throw std::runtime_error("SuperLU reuse run reported an invalid argument, info=" +
                                 std::to_string(reuseInfo));
    }
    if (reuseInfo > n + 1) {
        throw std::runtime_error("SuperLU reuse run failed with info=" + std::to_string(reuseInfo));
    }
    if (reuseInfo > 0 && reuseInfo <= n) {
        throw std::runtime_error("SuperLU reuse run detected a singular factor, info=" +
                                 std::to_string(reuseInfo));
    }

    result.times.fact_reuse = secondStat.utime[FACT];
    result.times.numeric = std::min(result.times.fact_reuse, result.times.fact_first);
    result.times.symbolic = result.times.etree +
                            std::max(0.0, result.times.fact_first - result.times.numeric);
    result.times.total = result.times.ordering + result.times.symbolic + result.times.numeric +
                         result.times.solve;

    return result;
}
