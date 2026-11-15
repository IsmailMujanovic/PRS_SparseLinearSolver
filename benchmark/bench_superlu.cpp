#include "../include/generators.hpp"
#include "../include/matrix.hpp"
#include "../include/md.hpp"
#include "../include/numeric.hpp"
#include "../include/permute.hpp"
#include "../include/symbolic.hpp"
#include "superlu_adapter.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct PhaseTimes {
    double ordering = 0.0;
    double symbolic = 0.0;
    double numeric = 0.0;
    double solve = 0.0;
    double total = 0.0;
};

struct OurResult {
    PhaseTimes times;
    double permute = 0.0;
    std::vector<Complex> solution;
    std::vector<std::size_t> permutation;
    std::size_t nnz_factor = 0;
    double factor_bytes = 0.0;
    double flops = 0.0;
};

struct Accuracy {
    double relative_residual = 0.0;
    double forward_error = 0.0;
};

struct MatrixCase {
    std::string name;
    Matrix matrix;
};

struct RandomCaseSpec {
    int n = 0;
    double density = 0.0;
    unsigned seed = 0;
};

struct Config {
    int repeat = 3;
    bool large = false;
    bool detailed = false;
    std::string csv;
    std::vector<int> grids;
    std::vector<RandomCaseSpec> randomCases;
};

struct SuperLURun {
    SuperLUResult result;
    Accuracy accuracy;
};

double elapsed(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

double median(std::vector<double> values) {
    if (values.empty()) {
        return 0.0;
    }
    std::sort(values.begin(), values.end());
    std::size_t middle = values.size() / 2;
    if (values.size() % 2 == 1) {
        return values[middle];
    }
    return 0.5 * (values[middle - 1] + values[middle]);
}

PhaseTimes medianPhases(const std::vector<PhaseTimes>& runs) {
    PhaseTimes result;
    std::vector<double> ordering;
    std::vector<double> symbolic;
    std::vector<double> numeric;
    std::vector<double> solve;
    std::vector<double> total;
    for (const auto& run : runs) {
        ordering.push_back(run.ordering);
        symbolic.push_back(run.symbolic);
        numeric.push_back(run.numeric);
        solve.push_back(run.solve);
        total.push_back(run.total);
    }
    result.ordering = median(ordering);
    result.symbolic = median(symbolic);
    result.numeric = median(numeric);
    result.solve = median(solve);
    result.total = median(total);
    return result;
}

std::vector<Complex> makeKnownSolution(std::size_t n) {
    std::vector<Complex> solution(n);
    for (std::size_t i = 0; i < n; i++) {
        solution[i] = Complex(static_cast<double>(i) + 1.0, 0.5 * static_cast<double>(i) - 1.0);
    }
    return solution;
}

std::vector<Complex> makeRightHandSide(const Matrix& A, const std::vector<Complex>& known) {
    std::vector<Complex> b(A.n, Complex(0.0, 0.0));
    for (std::size_t column = 0; column < A.n; column++) {
        for (std::size_t index = A.col_ptr[column]; index < A.col_ptr[column + 1]; index++) {
            b[A.row_indices[index]] += A.values[index] * known[column];
        }
    }
    return b;
}

double vectorInfinityNorm(const std::vector<Complex>& values) {
    double norm = 0.0;
    for (const auto& value : values) {
        norm = std::max(norm, std::abs(value));
    }
    return norm;
}

double matrixInfinityNorm(const Matrix& A) {
    std::vector<double> rowSums(A.n, 0.0);
    for (std::size_t column = 0; column < A.n; column++) {
        for (std::size_t index = A.col_ptr[column]; index < A.col_ptr[column + 1]; index++) {
            rowSums[A.row_indices[index]] += std::abs(A.values[index]);
        }
    }
    return *std::max_element(rowSums.begin(), rowSums.end());
}

Accuracy measureAccuracy(const Matrix& A, const std::vector<Complex>& b, const std::vector<Complex>& known,
                         const std::vector<Complex>& solution) {
    std::vector<Complex> residual = b;
    for (std::size_t column = 0; column < A.n; column++) {
        for (std::size_t index = A.col_ptr[column]; index < A.col_ptr[column + 1]; index++) {
            residual[A.row_indices[index]] -= A.values[index] * solution[column];
        }
    }

    double denominator = matrixInfinityNorm(A) * vectorInfinityNorm(solution) + vectorInfinityNorm(b);
    Accuracy accuracy;
    accuracy.relative_residual =
        denominator > 0.0 ? vectorInfinityNorm(residual) / denominator : vectorInfinityNorm(residual);

    double forward = 0.0;
    for (std::size_t i = 0; i < A.n; i++) {
        forward = std::max(forward, std::abs(solution[i] - known[i]));
    }
    accuracy.forward_error = forward;
    return accuracy;
}

double estimateNumericFlops(const SymbolicFactorization& symbolic) {
    std::size_t n = symbolic.n;
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> rowLists(n);
    for (std::size_t column = 0; column < n; column++) {
        for (std::size_t index = symbolic.col_ptr[column]; index < symbolic.col_ptr[column + 1]; index++) {
            rowLists[symbolic.row_indices[index]].push_back({column, index});
        }
    }

    std::size_t complexMultiplyAdds = 0;
    for (std::size_t column = 0; column < n; column++) {
        for (const auto& entry : rowLists[column]) {
            complexMultiplyAdds += symbolic.col_ptr[entry.first + 1] - entry.second - 1;
        }
    }

    return 8.0 * static_cast<double>(complexMultiplyAdds) +
           3.0 * static_cast<double>(symbolic.row_indices.size()) + 6.0 * static_cast<double>(n);
}

OurResult runOurSolver(const Matrix& A, const std::vector<Complex>& b, int repeat) {
    std::vector<PhaseTimes> runs;
    std::vector<double> permutes;
    OurResult last;

    for (int repetition = 0; repetition < repeat; repetition++) {
        auto start = std::chrono::steady_clock::now();
        auto permutation = MD(A);
        double ordering = elapsed(start);

        auto permutedStart = std::chrono::steady_clock::now();
        auto permuted = permute(A, permutation);
        double permutationTime = elapsed(permutedStart);

        auto symbolicStart = std::chrono::steady_clock::now();
        auto symbolic = symbolicFactorization(permuted);
        double symbolicTime = elapsed(symbolicStart);

        auto numericStart = std::chrono::steady_clock::now();
        auto factor = numericFactorization(permuted, symbolic);
        double numericTime = elapsed(numericStart);

        auto solveStart = std::chrono::steady_clock::now();
        auto solution = solve(factor, permutation, b);
        double solveTime = elapsed(solveStart);

        PhaseTimes phases;
        phases.ordering = ordering;
        phases.symbolic = permutationTime + symbolicTime;
        phases.numeric = numericTime;
        phases.solve = solveTime;
        phases.total = ordering + phases.symbolic + numericTime + solveTime;
        runs.push_back(phases);
        permutes.push_back(permutationTime);

        last.permutation = std::move(permutation);
        last.solution = std::move(solution);
        last.nnz_factor = factor.lower.values.size() + A.n;
        last.factor_bytes = static_cast<double>(factor.lower.values.size() * sizeof(Complex) +
                                                A.n * sizeof(double));
        last.flops = estimateNumericFlops(symbolic);
    }

    last.times = medianPhases(runs);
    last.permute = median(permutes);
    return last;
}

SuperLUPhaseTimes medianSuperLUPhases(const std::vector<SuperLUPhaseTimes>& runs) {
    SuperLUPhaseTimes result;
    std::vector<double> ordering;
    std::vector<double> symbolic;
    std::vector<double> numeric;
    std::vector<double> solve;
    std::vector<double> total;
    std::vector<double> etree;
    std::vector<double> factFirst;
    std::vector<double> factReuse;
    std::vector<double> wall;
    for (const auto& run : runs) {
        ordering.push_back(run.ordering);
        symbolic.push_back(run.symbolic);
        numeric.push_back(run.numeric);
        solve.push_back(run.solve);
        total.push_back(run.total);
        etree.push_back(run.etree);
        factFirst.push_back(run.fact_first);
        factReuse.push_back(run.fact_reuse);
        wall.push_back(run.wall);
    }
    result.ordering = median(ordering);
    result.symbolic = median(symbolic);
    result.numeric = median(numeric);
    result.solve = median(solve);
    result.total = median(total);
    result.etree = median(etree);
    result.fact_first = median(factFirst);
    result.fact_reuse = median(factReuse);
    result.wall = median(wall);
    return result;
}

SuperLURun runSuperLURepeated(const Matrix& A, const std::vector<Complex>& b, const std::vector<Complex>& known,
                              const SuperLUOptions& options, int repeat) {
    (void)runSuperLU(A, b, options);
    std::vector<SuperLUPhaseTimes> runs;
    SuperLURun last;
    for (int repetition = 0; repetition < repeat; repetition++) {
        last.result = runSuperLU(A, b, options);
        runs.push_back(last.result.times);
    }
    last.result.times = medianSuperLUPhases(runs);
    last.accuracy = measureAccuracy(A, b, known, last.result.solution);
    return last;
}

struct FillResult {
    std::size_t nnz_factor = 0;
    double seconds = 0.0;
};

FillResult fillUnderOrdering(const Matrix& A, const std::vector<std::size_t>& ordering) {
    auto start = std::chrono::steady_clock::now();
    auto permuted = permute(A, ordering);
    auto symbolic = symbolicFactorization(permuted);
    FillResult result;
    result.seconds = elapsed(start);
    result.nnz_factor = symbolic.row_indices.size() + A.n;
    return result;
}

std::vector<MatrixCase> buildCases(const Config& config) {
    std::vector<MatrixCase> cases;
    cases.push_back({"tridiagonal-300", tridiagonalHPD(300)});
    cases.push_back({"random-60", randomHPD(60, 0.08, 1)});
    cases.push_back({"random-300", randomHPD(300, 0.02, 2)});
    cases.push_back({"laplacian2d-20", complexLaplacian2D(20)});
    cases.push_back({"laplacian2d-30", complexLaplacian2D(30)});
    if (config.large) {
        cases.push_back({"random-900", randomHPD(900, 2.0 / 900.0, 3)});
        cases.push_back({"laplacian2d-40", complexLaplacian2D(40)});
    }
    for (int grid : config.grids) {
        if (grid < 1) {
            throw std::runtime_error("grid size must be positive");
        }
        cases.push_back({"custom-grid-" + std::to_string(grid),
                         complexLaplacian2D(static_cast<std::size_t>(grid))});
    }
    for (const auto& spec : config.randomCases) {
        if (spec.n < 1 || spec.density <= 0.0 || spec.density > 1.0) {
            throw std::runtime_error("random matrix requires n >= 1 and density in (0, 1]");
        }
        cases.push_back({"custom-random-" + std::to_string(spec.n),
                         randomHPD(static_cast<std::size_t>(spec.n), spec.density, spec.seed)});
    }
    return cases;
}

std::string formatRatio(double numerator, double denominator) {
    if (numerator == 0.0 && denominator == 0.0) {
        return "1.000";
    }
    if (denominator == 0.0) {
        return "inf";
    }
    double value = numerator / denominator;
    if (!std::isfinite(value)) {
        return "inf";
    }
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(3) << value;
    return stream.str();
}

std::string formatNumber(double value, int precision = 3) {
    std::ostringstream stream;
    if (!std::isfinite(value)) {
        return "inf";
    }
    stream << std::fixed << std::setprecision(precision) << value;
    return stream.str();
}

std::string scientific(double value) {
    std::ostringstream stream;
    stream << std::scientific << std::setprecision(2) << value;
    return stream.str();
}

void printSimpleTable(const std::string& name, const Matrix& A, const OurResult& ours,
                      const SuperLURun& reference) {
    std::cout << "\n" << name << "  n=" << A.n << "  nnz(A)=" << A.values.size() << "\n";
    std::cout << std::left << std::setw(24) << "Phase" << std::right << std::setw(14) << "ours (ms)"
              << std::setw(16) << "SuperLU (ms)" << std::setw(22) << "ours - SuperLU (ms)"
              << std::setw(20) << "ours / SuperLU (x)" << "\n";

    auto printRow = [](const std::string& label, double ourSeconds, double superluSeconds) {
        std::cout << std::left << std::setw(24) << label << std::right << std::setw(14)
                  << formatNumber(1000.0 * ourSeconds) << std::setw(16)
                  << formatNumber(1000.0 * superluSeconds) << std::setw(22)
                  << formatNumber(1000.0 * (ourSeconds - superluSeconds)) << std::setw(20)
                  << formatRatio(ourSeconds, superluSeconds) << "\n";
    };

    printRow("1 Ordering", ours.times.ordering, reference.result.times.ordering);
    printRow("2 Symbolic factorization", ours.times.symbolic, reference.result.times.symbolic);
    printRow("3 Numeric factorization", ours.times.numeric, reference.result.times.numeric);
    printRow("4 Triangular solution", ours.times.solve, reference.result.times.solve);
    printRow("Total", ours.times.total, reference.result.times.total);
}

void printPhaseTable(const std::string& name, const Matrix& A, const OurResult& our, const SuperLURun& reference) {
    const SuperLUResult& superlu = reference.result;
    std::cout << "\nMatrix: " << name << "  n=" << A.n << "  nnz(A)=" << A.values.size() << "\n";
    std::cout << std::left << std::setw(24) << "Phase" << std::right << std::setw(12) << "ours ms"
              << std::setw(12) << "SLU ms" << std::setw(11) << "slowdown" << std::setw(10) << "our %"
              << std::setw(11) << "gap ms" << std::setw(9) << "gap %" << "\n";

    auto printRow = [&](const std::string& label, double ourTime, double superluTime, double totalGap) {
        double gap = 1000.0 * (ourTime - superluTime);
        double ourShare = our.times.total > 0.0 ? 100.0 * ourTime / our.times.total : 0.0;
        double gapShare = totalGap != 0.0 ? 100.0 * (ourTime - superluTime) / totalGap : 0.0;
        std::cout << std::left << std::setw(24) << label << std::right << std::setw(12)
                  << formatNumber(1000.0 * ourTime) << std::setw(12) << formatNumber(1000.0 * superluTime)
                  << std::setw(11) << formatRatio(ourTime, superluTime) << std::setw(10)
                  << formatNumber(ourShare, 1) << std::setw(11) << formatNumber(gap)
                  << std::setw(9) << formatNumber(gapShare, 1) << "\n";
    };

    double totalGap = our.times.total - superlu.times.total;
    printRow("1 Ordering", our.times.ordering, superlu.times.ordering, totalGap);
    printRow("2 Symbolic (+permute)", our.times.symbolic, superlu.times.symbolic, totalGap);
    printRow("3 Numeric", our.times.numeric, superlu.times.numeric, totalGap);
    printRow("4 Triangular solve", our.times.solve, superlu.times.solve, totalGap);
    printRow("Total", our.times.total, superlu.times.total, totalGap);

    std::cout << "  ours: permute=" << formatNumber(1000.0 * our.permute) << " ms, nnz(L)=" << our.nnz_factor
              << ", factor=" << formatNumber(our.factor_bytes / (1024.0 * 1024.0)) << " MiB, flops="
              << formatNumber(our.flops, 0) << "\n";
    std::cout << "  superlu: nnz(L+U-n)=" << superlu.nnz_factor << ", chol-equivalent="
              << (superlu.nnz_l + superlu.nnz_u) / 2 << ", factor="
              << formatNumber(superlu.factor_bytes / (1024.0 * 1024.0)) << " MiB, flops="
              << formatNumber(superlu.flops, 0) << ", info=" << superlu.info
              << ", split=" << (superlu.split_valid ? "estimated" : "raw") << "\n";
    double ourRate = our.times.numeric > 0.0 ? our.flops / (1e9 * our.times.numeric) : 0.0;
    double superluRate =
        superlu.times.numeric > 0.0 ? superlu.flops / (1e9 * superlu.times.numeric) : 0.0;
    std::cout << "  numeric rate: ours=" << formatNumber(ourRate) << " GFLOP/s, SuperLU="
              << formatNumber(superluRate) << " GFLOP/s\n";
}

void printModeTable(const std::string& name, const OurResult& ours, const Accuracy& ourAccuracy,
                    std::size_t nnzA, const std::vector<std::pair<std::string, SuperLURun>>& modes) {
    std::cout << "\nSuperLU variants for " << name << "\n";
    std::cout << std::left << std::setw(14) << "solver/mode" << std::right << std::setw(10) << "order"
              << std::setw(11) << "symb est" << std::setw(11) << "numeric" << std::setw(10) << "solve"
              << std::setw(10) << "total" << std::setw(11) << "nnz(LU)" << std::setw(9) << "fill"
              << std::setw(12) << "relres" << std::setw(8) << "info" << "\n";

    auto printRow = [&](const std::string& label, const auto& times, std::size_t nnzFactor,
                        const Accuracy& accuracy, int info) {
        double fill = static_cast<double>(nnzFactor) / static_cast<double>(nnzA);
        std::cout << std::left << std::setw(14) << label << std::right << std::setw(10)
                  << formatNumber(1000.0 * times.ordering) << std::setw(11) << formatNumber(1000.0 * times.symbolic)
                  << std::setw(11) << formatNumber(1000.0 * times.numeric) << std::setw(10)
                  << formatNumber(1000.0 * times.solve) << std::setw(10) << formatNumber(1000.0 * times.total)
                  << std::setw(11) << nnzFactor << std::setw(9) << formatNumber(fill, 2)
                  << std::setw(12) << scientific(accuracy.relative_residual) << std::setw(8) << info
                  << "\n";
    };

    printRow("ours/MD", ours.times, ours.nnz_factor, ourAccuracy, 0);
    for (const auto& entry : modes) {
        printRow(entry.first, entry.second.result.times, entry.second.result.nnz_factor,
                 entry.second.accuracy, entry.second.result.info);
    }
}

void printOrderingTable(const MatrixCase& testCase, const OurResult& ours,
                        const std::vector<std::pair<std::string, SuperLUOrderingResult>>& orderings) {
    std::cout << "\nOrdering quality for " << testCase.name << " (same symbolic code for every ordering)\n";
    std::cout << std::left << std::setw(16) << "ordering" << std::right << std::setw(12) << "time ms"
              << std::setw(12) << "nnz(L)" << std::setw(12) << "fill ratio" << std::setw(12)
              << "symb ms" << "\n";

    std::cout << std::left << std::setw(16) << "our MD" << std::right << std::setw(12)
              << formatNumber(1000.0 * ours.times.ordering) << std::setw(12) << ours.nnz_factor << std::setw(12)
              << formatNumber(static_cast<double>(ours.nnz_factor) /
                              static_cast<double>(testCase.matrix.values.size()))
              << std::setw(12) << formatNumber(1000.0 * ours.times.symbolic) << "\n";

    for (const auto& entry : orderings) {
        FillResult fill = fillUnderOrdering(testCase.matrix, entry.second.permutation);
        std::cout << std::left << std::setw(16) << entry.first << std::right << std::setw(12)
                  << formatNumber(entry.second.seconds * 1000.0) << std::setw(12) << fill.nnz_factor
                  << std::setw(12)
                  << formatNumber(static_cast<double>(fill.nnz_factor) /
                                  static_cast<double>(testCase.matrix.values.size()))
                  << std::setw(12) << formatNumber(fill.seconds * 1000.0) << "\n";
    }
}

void printDiagnosis(const MatrixCase& testCase, const OurResult& ours, const SuperLURun& reference,
                    const SuperLUOrderingResult& mmdOrdering, const SuperLUOrderingResult& colamdOrdering,
                    const SuperLURun& sameOrder) {
    double gaps[4] = {
        ours.times.ordering - reference.result.times.ordering,
        ours.times.symbolic - reference.result.times.symbolic,
        ours.times.numeric - reference.result.times.numeric,
        ours.times.solve - reference.result.times.solve};
    const char* labels[4] = {"ordering", "symbolic", "numeric", "triangular solve"};
    int largestGap = 0;
    for (int i = 1; i < 4; i++) {
        if (gaps[i] > gaps[largestGap]) {
            largestGap = i;
        }
    }

    double shares[4] = {ours.times.ordering, ours.times.symbolic, ours.times.numeric, ours.times.solve};
    int largestShare = 0;
    for (int i = 1; i < 4; i++) {
        if (shares[i] > shares[largestShare]) {
            largestShare = i;
        }
    }

    double mdFill = static_cast<double>(ours.nnz_factor);
    FillResult mmdFill = fillUnderOrdering(testCase.matrix, mmdOrdering.permutation);
    FillResult colamdFill = fillUnderOrdering(testCase.matrix, colamdOrdering.permutation);

    std::cout << "Diagnosis for " << testCase.name << ":\n";
    std::cout << "  largest phase gap: " << labels[largestGap] << " ("
              << formatNumber(1000.0 * gaps[largestGap]) << " ms)\n";
    std::cout << "  largest share of our total: " << labels[largestShare] << " ("
              << formatNumber(100.0 * shares[largestShare] / ours.times.total, 1) << "%)\n";
    std::cout << "  fill with our MD: " << mdFill << ", with SuperLU MMD: " << mmdFill.nnz_factor
              << ", with SuperLU COLAMD: " << colamdFill.nnz_factor << "\n";
    std::cout << "  same-ordering numeric: ours=" << formatNumber(1000.0 * ours.times.numeric)
              << " ms, SuperLU=" << formatNumber(1000.0 * sameOrder.result.times.numeric)
              << " ms, slowdown=" << formatRatio(ours.times.numeric, sameOrder.result.times.numeric)
              << "\n";
}

void writeSimpleCsvRows(std::ofstream& csv, const std::string& matrix, const Matrix& A,
                        const OurResult& ours, const SuperLUResult& superlu) {
    auto writeRow = [&](const std::string& phase, double ourSeconds, double superluSeconds) {
        double difference = (ourSeconds - superluSeconds) * 1000.0;
        double ratio = superluSeconds > 0.0 ? ourSeconds / superluSeconds : 0.0;
        csv << matrix << ',' << A.n << ',' << phase << ',' << ourSeconds * 1000.0 << ','
            << superluSeconds * 1000.0 << ',' << difference << ',' << ratio << '\n';
    };

    writeRow("1 Ordering", ours.times.ordering, superlu.times.ordering);
    writeRow("2 Symbolic factorization", ours.times.symbolic, superlu.times.symbolic);
    writeRow("3 Numeric factorization", ours.times.numeric, superlu.times.numeric);
    writeRow("4 Triangular solution", ours.times.solve, superlu.times.solve);
    writeRow("Total", ours.times.total, superlu.times.total);
}

Config parseArguments(int argc, char** argv) {
    Config config;
    for (int i = 1; i < argc; i++) {
        std::string argument = argv[i];
        if (argument == "--repeat") {
            if (i + 1 >= argc) {
                throw std::runtime_error("--repeat requires a value");
            }
            config.repeat = std::stoi(argv[++i]);
            if (config.repeat < 1) {
                throw std::runtime_error("--repeat must be at least 1");
            }
        } else if (argument == "--large") {
            config.large = true;
        } else if (argument == "--detailed") {
            config.detailed = true;
        } else if (argument == "--simple") {
            config.detailed = false;
        } else if (argument == "--csv") {
            if (i + 1 >= argc) {
                throw std::runtime_error("--csv requires a path");
            }
            config.csv = argv[++i];
        } else if (argument == "--grid") {
            if (i + 1 >= argc) {
                throw std::runtime_error("--grid requires a value");
            }
            config.grids.push_back(std::stoi(argv[++i]));
        } else if (argument == "--random") {
            if (i + 1 >= argc) {
                throw std::runtime_error("--random requires N:DENSITY:SEED");
            }
            std::string specification = argv[++i];
            auto first = specification.find(':');
            auto second = specification.find(':', first == std::string::npos ? 0 : first + 1);
            if (first == std::string::npos || second == std::string::npos) {
                throw std::runtime_error("--random requires N:DENSITY:SEED");
            }
            RandomCaseSpec randomCase;
            randomCase.n = std::stoi(specification.substr(0, first));
            randomCase.density = std::stod(specification.substr(first + 1, second - first - 1));
            randomCase.seed = static_cast<unsigned>(std::stoul(specification.substr(second + 1)));
            config.randomCases.push_back(randomCase);
        } else if (argument == "--help") {
            std::cout << "usage: bench_superlu [--repeat N] [--large] [--detailed] [--csv FILE] "
                         "[--grid K] [--random N:DENSITY:SEED]\n";
            std::exit(0);
        } else {
            throw std::runtime_error("unknown argument: " + argument);
        }
    }
    return config;
}

}

int main(int argc, char** argv) {
    try {
        Config config = parseArguments(argc, argv);
        auto cases = buildCases(config);

        std::ofstream csv;
        if (!config.csv.empty()) {
            csv.open(config.csv);
            if (!csv) {
                throw std::runtime_error("cannot open CSV output: " + config.csv);
            }
            csv << "matrix,n,phase,ours_ms,superlu_ms,ours_minus_superlu_ms,ours_over_superlu_x\n";
        }

        bool ok = true;
        for (const auto& testCase : cases) {
            try {
                const Matrix& A = testCase.matrix;
                auto known = makeKnownSolution(A.n);
                auto b = makeRightHandSide(A, known);

                OurResult ours = runOurSolver(A, b, config.repeat);
                Accuracy ourAccuracy = measureAccuracy(A, b, known, ours.solution);

                SuperLUOptions hpdOptions;
                hpdOptions.mode = SuperLUMode::Hpd;
                SuperLURun hpd = runSuperLURepeated(A, b, known, hpdOptions, config.repeat);

                if (csv) {
                    writeSimpleCsvRows(csv, testCase.name, A, ours, hpd.result);
                }

                if (config.detailed) {
                    SuperLUOptions defaultOptions;
                    defaultOptions.mode = SuperLUMode::Default;
                    SuperLURun standard = runSuperLURepeated(A, b, known, defaultOptions, config.repeat);

                    SuperLUOptions sameOrderOptions;
                    sameOrderOptions.mode = SuperLUMode::SameOrder;
                    sameOrderOptions.ordering = &ours.permutation;
                    SuperLURun sameOrder = runSuperLURepeated(A, b, known, sameOrderOptions, config.repeat);

                    auto mmdOrdering = superluOrdering(A, SuperLUOrderingMethod::MmdAtPlusA);
                    auto colamdOrdering = superluOrdering(A, SuperLUOrderingMethod::Colamd);

                    printPhaseTable(testCase.name, A, ours, hpd);
                    printModeTable(testCase.name, ours, ourAccuracy, A.values.size(),
                                   {{"default", standard}, {"hpd", hpd}, {"same-order", sameOrder}});
                    printOrderingTable(testCase, ours,
                                       {{"SuperLU MMD", mmdOrdering}, {"SuperLU COLAMD", colamdOrdering}});
                    printDiagnosis(testCase, ours, hpd, mmdOrdering, colamdOrdering, sameOrder);

                    std::cout << "  accuracy: ours relres=" << std::scientific << std::setprecision(3)
                              << ourAccuracy.relative_residual << " fwd=" << ourAccuracy.forward_error
                              << ", SLU relres=" << hpd.accuracy.relative_residual
                              << " fwd=" << hpd.accuracy.forward_error << std::defaultfloat << "\n";
                } else {
                    printSimpleTable(testCase.name, A, ours, hpd);
                }

                if (ourAccuracy.relative_residual > 1e-8 || hpd.accuracy.relative_residual > 1e-8) {
                    ok = false;
                }
            } catch (const std::exception& error) {
                std::cout << "ERROR in " << testCase.name << ": " << error.what() << "\n";
                ok = false;
            }
            std::cout << std::flush;
        }

        return ok ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "benchmark failed: " << error.what() << "\n";
        return 1;
    }
}
