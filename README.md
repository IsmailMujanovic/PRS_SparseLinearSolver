# PRS Sparse Linear Solver

Ovaj projekat je početna implementacija C++ rješenja za rad sa rijetkim matricama i redukciju grafova koristeći Minimum Degree (MD) algoritam.

## Šta je urađeno do sada

U trenutnoj verziji projekat sadrži:

- definiciju strukture `Matrix` za čuvanje rijetke matrice,
- podatke o vrijednostima matrice, indeksima redova i pointerima kolona,
- implementaciju `MD_Graph` i funkcije `MD(...)` za traženje permutacije čvorova,
- jednostavan primjer u `main.cpp` koji kreira malu matricu i ispisuje rezultat.

Ovo je još osnovna verzija, a cilj je da se kasnije proširi za efikasnije rješavanje linearnih sistema sa rijetkim matricama.

## Struktura projekta

- `include/matrix.hpp` - deklaracija strukture matrice
- `include/md.hpp` - deklaracija grafa i MD algoritma
- `src/matrix.cpp` - konstruktor matrice
- `src/md.cpp` - implementacija Minimum Degree algoritma
- `include/permute.hpp`, `src/permute.cpp` - permutacija matrice (donji trougao od P A P^H)
- `include/symbolic.hpp`, `src/symbolic.cpp` - simbolička faktorizacija (struktura L i elimination tree)
- `include/numeric.hpp`, `src/numeric.cpp` - numerička Cholesky faktorizacija (A = L L^H) i rješavanje sistema
- `src/main.cpp` - primjer korištenja i provjera tačnosti

## Kako pokrenuti projekat

### Windows (PowerShell)

U direktoriju projekta izvršite:

```powershell
g++ -std=c++17 src/main.cpp src/matrix.cpp src/md.cpp src/permute.cpp src/symbolic.cpp src/numeric.cpp -Iinclude -o program.exe
./program.exe
```

Ako je program već kompajliran, dovoljno je pokrenuti:

```powershell
./program.exe
```

### Linux / macOS

U direktoriju projekta izvršite:

```bash
g++ -std=c++17 src/main.cpp src/matrix.cpp src/md.cpp src/permute.cpp src/symbolic.cpp src/numeric.cpp -Iinclude -o program
./program
```

Ako je program već kompajliran, dovoljno je pokrenuti:

```bash
./program
```

> Ako se pojavi greška da g++ nije instaliran, instalirajte C++ kompajler za vaš operativni sistem (npr. `build-essential` na Ubuntu/Debian, ili Xcode Command Line Tools na macOS).

## Šta program radi

Program za nekoliko kompleksnih Hermitian pozitivno-definitnih matrica (mala 4x4 i dvije slučajne rijetke, n=60 i n=300) izvršava cijeli tok:

1. `MD(A)` - permutacija (Minimum Degree),
2. `permute(A, P)` - donji trougao od P A P^H,
3. `symbolicFactorization(...)` - struktura L i elimination tree (marker niz, O(eta(L))),
4. `numericFactorization(...)` - left-looking column Cholesky (cmod/cdiv) uz kompleksnu konjugaciju,
5. `solve(...)` - rješavanje A x = b.

Za svaku matricu ispisuje `|PAP^H - L L^H|` i grešku rješenja, te na kraju `ALL OK` ako su obje greške manje od 1e-10.
