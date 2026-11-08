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
- `src/main.cpp` - primjer korištenja

## Kako pokrenuti projekat

### Windows (PowerShell)

U direktoriju projekta izvršite:

```powershell
g++ src/main.cpp src/matrix.cpp src/md.cpp -Iinclude -o program.exe
./program.exe
```

Ako je program već kompajliran, dovoljno je pokrenuti:

```powershell
./program.exe
```

### Linux / macOS

U direktoriju projekta izvršite:

```bash
g++ src/main.cpp src/matrix.cpp src/md.cpp -Iinclude -o program
./program
```

Ako je program već kompajliran, dovoljno je pokrenuti:

```bash
./program
```

> Ako se pojavi greška da g++ nije instaliran, instalirajte C++ kompajler za vaš operativni sistem (npr. `build-essential` na Ubuntu/Debian, ili Xcode Command Line Tools na macOS).

## Šta program radi

Program trenutno kreira malu test matricu i poziva funkciju `MD(matrix)`. Rezultat je lista permutacije čvorova koja se ispisuje na konzolu.

To je početni korak ka razvoju algoritama za redukciju grafa i kasnije efikasnijeg rješavanja sistema linearnih jednadžbi.
