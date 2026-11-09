# Sparse Cholesky pipeline: what we did, step by step

Goal: solve `A x = b` for a sparse complex Hermitian positive-definite (HPD) matrix `A`, following the Heath, Ng and Peyton paper (Sections 2.2–2.4).

The pipeline is:

```
A --MD--> P --permute--> PAP^H --symbolic--> structure of L --numeric--> L --solve--> x
```

Status: the code compiles and has been verified end-to-end. The built-in tests print `ALL OK`; see §8 for what is checked and §9 for the build command.

---

## 1. Starting point

The repo only had:

- `Matrix` (`include/matrix.hpp`): CSC storage with `values` (complex), `row_indices`, `col_ptr`.
- `MD(...)` (`src/md.cpp`): Minimum Degree ordering that returns a permutation `P`, where `P[i]` is the original vertex eliminated i-th.

Missing: permuting `A`, symbolic factorization, numeric factorization and the triangular solve.

## 2. Review of the first proposed code

The first draft followed the paper's idea, but had problems that we fixed before implementing:

| Problem | Why it is wrong | Fix |
|---|---|---|
| `std::set` for column structures | Costs O(log n) per insert, so overall O(eta(L) log n) | Marker array `mark[r] == j` plus sorting only at the end of each column |
| `permute` summed duplicate entries | `Matrix` stores both triangles. `(r,c)` and its mirror `(c,r)` land on the same lower position with the same value, so the sum doubles every off-diagonal | Read only the lower triangle of `A` (`r >= c`); each entry is placed exactly once |
| Numeric step lost the diagonal of `L` | The symbolic structure excludes the diagonal and `diag_sqrt` was never stored | `NumericFactor` keeps a separate `diag` vector |
| `for k in 0..j` with `L[k].count(j)` | Scans all earlier columns for every column, which is O(n^2) or worse | `rowlist[j]` lists only the columns that have a nonzero in row `j` |
| Column rows unsorted, `if (r > j)` tests | Each `l_jk` update rescanned the whole column k | Sort rows once; start the inner loop right after `l_jk` |
| `MD_Graph` skipped the first entry of each column | Assumes it is the diagonal. With both triangles stored this is only true for column 0, so it created self-loops and dropped edges | Skip entries with `row == i` instead |
| Pasted code redefined `Matrix`, `Complex`, `MD_Graph` | Clashes with the repo headers | Split into separate headers and sources under `include/` and `src/` |

## 3. Step 0: fix the MD graph (`src/md.cpp`)

In the `MD_Graph` constructor the loop used to start at `col_ptr[i] + 1`. It now starts at `col_ptr[i]` and skips entries where `row_indices[j] == i`. This gives correct adjacency lists, so the permutation feeding everything downstream is correct. `eliminateVertex` was left unchanged.

## 4. Step 1: permute (`include/permute.hpp`, `src/permute.cpp`)

`Matrix permute(const Matrix& A, const std::vector<std::size_t>& P)` builds the lower triangle (diagonal included) of `P A P^H`.

1. Build the inverse permutation: `invP[P[i]] = i`.
2. Loop over the entries of `A`, skipping every `r < c`. The upper triangle is the conjugate mirror of the lower one, so it carries no new information.
3. Map `(r, c)` to `(newR, newC) = (invP[r], invP[c])`.
   - If `newR >= newC`, store the value in column `newC` at row `newR`.
   - Otherwise the entry landed in the upper triangle. Swap the indices and store `conj(value)`, because for a Hermitian matrix `a_cr = conj(a_rc)`.
4. Sort each column by row and write it out as a CSC `Matrix`.

## 5. Step 2: symbolic factorization (`include/symbolic.hpp`, `src/symbolic.cpp`)

Computes the sparsity pattern of `L` and the elimination tree without touching any numeric values (paper Section 2.2).

```cpp
struct SymbolicFactorization {
    std::size_t n;
    std::vector<std::size_t> col_ptr;
    std::vector<std::size_t> row_indices;
    std::vector<std::size_t> parent;
};
```

`col_ptr` and `row_indices` store the structure of `L` strictly below the diagonal, with rows sorted. `parent[j] == n` marks a root of the elimination tree.

For every column `j`, from 0 to n-1:

1. Add every row `r > j` that appears in column `j` of `PAPt`.
2. Add every row `r > j` from the columns of all children `c` of `j` (columns whose parent is `j`). Fill-in from a child spills into its parent.
3. The smallest row in the column is the parent `p(j)`. Record it and add `j` to `children[p(j)]`.
4. Sort the rows of column `j`.

A marker array makes the "add" step O(1) with no duplicates: row `r` is added only if `mark[r] != j`, then `mark[r] = j`. Total cost is O(eta(A) + eta(L)), plus the per-column sort.

## 6. Step 3: numeric factorization (`include/numeric.hpp`, `src/numeric.cpp`)

Left-looking column Cholesky `PAP^H = L L^H` (paper Section 2.3).

```cpp
struct NumericFactor {
    std::vector<double> diag;
    Matrix lower;
};
```

`diag` is the real diagonal of `L`. `lower` is `L` strictly below the diagonal, reusing the symbolic pattern, so no fill-in slot has to be allocated later.

Setup: build `rowlist[r]` once. It contains `(k, pos)` for every column `k < r` with `l_rk != 0`, where `pos` is the position of `l_rk` in `lower.values`.

For every column `j`:

1. **Scatter** column `j` of `PAPt` into a dense work vector `w` of size n.
2. **cmod**: for every `(k, pos)` in `rowlist[j]`:
   - `w[j] -= |l_jk|^2`
   - for every later entry `idx` in column `k`: `w[row] -= l_{row,k} * conj(l_jk)`. Rows are sorted, so all of them are greater than `j`.

   The `conj` is what makes it work for Hermitian matrices (`A = L L^H`).
3. **cdiv**: the pivot is `w[j].real()`. If it is `<= 0` the matrix is not positive-definite, so throw `std::runtime_error`. Otherwise `diag[j] = sqrt(pivot)`.
4. **Gather**: divide each `w[row]` of the column pattern by `diag[j]`, store it in `lower.values`, and reset `w` to zero for the next column.

Because we only touch entries that are actually nonzero in `L`, the work is proportional to the number of operations Cholesky really needs.

## 7. Step 4: triangular solution (`solve` in `src/numeric.cpp`)

Paper Section 2.4. Once `PAP^H = L L^H` is available, solving `A x = b` reduces to two triangular systems on the permuted right-hand side:

```
L y = P b,      L^H z = y,      x = P^H z.
```

`std::vector<Complex> solve(const NumericFactor& L, const std::vector<std::size_t>& P, const std::vector<Complex>& b)`

1. Permute the right-hand side: `y[i] = b[P[i]]`.
2. Forward substitution `L y = P b`, column by column: divide by `diag[j]`, then subtract the column from the rows below.
3. Backward substitution `L^H z = y`, columns from last to first: subtract `conj(l) * y[row]` for each entry of the column, then divide by `diag[j]`.
4. Undo the permutation: `x[P[i]] = y[i]`.

Both substitutions walk the CSC arrays of `L` directly, so only nonzeros of `L` are touched. The column order respects the elimination tree `T(A)` implicitly: Cholesky guarantees `parent[j] > j`, so ascending order for `L` and descending order for `L^H` are valid topological orders of the tree, and no explicit tree traversal is needed in the serial algorithm. The work is proportional to `eta(L)` nonzeros rather than the dense `n(n-1)/2`, which is why the triangular solution is usually a small fraction of the total serial time (paper Section 2.4). Parallel triangular solution (paper Section 3.5) is not implemented.

## 8. Verification (`src/main.cpp`)

`check(A, name)` runs the whole pipeline and tests two things:

1. **Factorization error**: build dense `PAP^H` and dense `L` (diagonal plus lower part), and compute the largest entry of `|PAP^H - L L^H|`.
2. **Solve error**: pick a known `x`, set `b = A x`, solve, and compute the largest entry of `|x_solved - x_known|`.

Both must be below `1e-10`. Test inputs:

- the 4x4 complex Hermitian tridiagonal matrix (both triangles stored),
- two random sparse Hermitian diagonally dominant matrices (diagonally dominant implies positive-definite) with n = 60 and n = 300, produced by `randomHPD`.

The program prints `ALL OK` and returns 0 if every check passes, otherwise `FAILED` and returns 1.

## 9. Build and run

```
g++ -std=c++17 -Wall -Wextra src/main.cpp src/matrix.cpp src/md.cpp src/permute.cpp src/symbolic.cpp src/numeric.cpp -Iinclude -o program.exe
./program.exe
```

Expected last line: `ALL OK`. If compilation or a check fails, the output will show which step is wrong (a large `|PAP^H-LL^H|` points at permute, symbolic or numeric; a large solve error with a good factorization points at `solve`).

## 10. File map

| File | Role |
|---|---|
| `include/matrix.hpp`, `src/matrix.cpp` | CSC `Matrix` (existing) |
| `include/md.hpp`, `src/md.cpp` | Minimum Degree ordering (diagonal bug fixed) |
| `include/permute.hpp`, `src/permute.cpp` | Lower triangle of `P A P^H` |
| `include/symbolic.hpp`, `src/symbolic.cpp` | Structure of `L` and elimination tree |
| `include/numeric.hpp`, `src/numeric.cpp` | Numeric Cholesky and triangular solution (`solve`) |
| `src/main.cpp` | Test inputs and verification |
| `README.md` | Short project overview and build commands |
