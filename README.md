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
- `include/validate.hpp`, `src/validate.cpp` - provjera da je ulazna matrica Hermitijska i pozitivno definitna
- `include/permute.hpp`, `src/permute.cpp` - permutacija matrice (donji trougao od P A P^H)
- `include/symbolic.hpp`, `src/symbolic.cpp` - simbolička faktorizacija (struktura L i elimination tree)
- `include/numeric.hpp`, `src/numeric.cpp` - numerička Cholesky faktorizacija (A = L L^H) i rješavanje sistema
- `include/generators.hpp` - generatori testnih matrica (slučajna HPD, tridijagonalna, kompleksni 2D Laplacian)
- `benchmark/` - SuperLU benchmark: adapter, četiri faze, CSV i analiza
- `run-benchmark.sh`, `run-benchmark.ps1`, `run-benchmark.bat` - pokretači benchmarka za macOS/Linux i Windows
- `CMakeLists.txt`, `Makefile` - build sistem za Windows (CMake) i macOS/Linux (Make)
- `docs/running.md` - uputstvo za pokretanje na svakoj platformi
- `src/main.cpp` - primjer korištenja i provjera tačnosti

## Kako pokrenuti projekat

### Windows (PowerShell)

U direktoriju projekta izvršite:

```powershell
g++ -std=c++17 src/main.cpp src/matrix.cpp src/md.cpp src/permute.cpp src/symbolic.cpp src/numeric.cpp src/validate.cpp -Iinclude -o program.exe
./program.exe
```

Ako je program već kompajliran, dovoljno je pokrenuti:

```powershell
./program.exe
```

### Linux / macOS

U direktoriju projekta izvršite:

```bash
g++ -std=c++17 src/main.cpp src/matrix.cpp src/md.cpp src/permute.cpp src/symbolic.cpp src/numeric.cpp src/validate.cpp -Iinclude -o program
./program
```

Ako je program već kompajliran, dovoljno je pokrenuti:

```bash
./program
```

> Ako se pojavi greška da g++ nije instaliran, instalirajte C++ kompajler za vaš operativni sistem (npr. `build-essential` na Ubuntu/Debian, ili Xcode Command Line Tools na macOS).

## Provjera ulazne matrice

Prije obrade svake matrice poziva se `validateHPD(A)` iz `include/validate.hpp`. Provjera ima dva dijela:

1. **Hermitijska simetrija** (`A = A^H`): za kompleksne matrice to znači `a[i][j] = conj(a[j][i])`. Dijagonala mora biti realna, a svaki upisani element mora se poklapati sa svojim konjugovanim parom. Elementi koji nisu upisani tretiraju se kao nula, pa provjera radi i kada je sačuvan samo donji trougao.
2. **Pozitivna definitnost** (`x^H A x > 0` za svako `x != 0`): provjerava se nezavisnom Cholesky faktorizacijom `A = L L^H`. Ako je matrica Hermitijska, ona je pozitivno definitna ako i samo ako svi pivoti te faktorizacije ostanu konačni i pozitivni.

Provjere koriste relativnu toleranciju (podrazumijevano `1e-12`), pa ne zavise od reda veličine elemenata. Ako matrica ne prođe provjeru, program ispisuje razlog (npr. koji element kvari simetriju ili na kojem indeksu je pivot nepravilan), preskače obradu te matrice i na kraju prijavljuje `FAILED`. Provjera je namjerno implementirana odvojeno od `numericFactorization`, kako greška u samom rješavaču ne bi mogla sakriti neispravan ulaz.

## Šta program radi

Program za nekoliko kompleksnih Hermitijskih pozitivno-definitnih matrica (mala 4x4 i dvije slučajne rijetke, n=60 i n=300) izvršava cijeli tok:

0. `validateHPD(A)` - provjera ulazne matrice (detalji u sekciji "Provjera ulazne matrice"),
1. `MD(A)` - permutacija (Minimum Degree),
2. `permute(A, P)` - donji trougao od P A P^H,
3. `symbolicFactorization(...)` - struktura L i elimination tree (marker niz, O(eta(L))),
4. `numericFactorization(...)` - left-looking column Cholesky (cmod/cdiv) uz kompleksnu konjugaciju,
5. `solve(...)` - rješavanje A x = b.

Za svaku ispravnu matricu ispisuje `|PAP^H - L L^H|` i grešku rješenja, te na kraju `ALL OK` ako su obje greške manje od 1e-10. Program dodatno testira odbijanje ne-Hermitijske matrice, Hermitijske nepozitivno-definitne matrice i matrice sa nerealnom dijagonalom.

## Benchmark sa SuperLU

Projekat sadrži i benchmark koji poredi sve četiri faze rješavača (ordering, simbolička faktorizacija, numerička faktorizacija i trougaono rješavanje) sa SuperLU 7.0.1. SuperLU se koristi kao referentna, visoko optimizovana implementacija; nije matematički optimalan rješavač i računa LU faktorizaciju, pa se poređenje popunjenosti i performansi prilagođava toj razlici.

Instalacija SuperLU na macOS-u:

```bash
brew install superlu
```

Najjednostavnije pokretanje na macOS-u i Linuxu:

```bash
./run-benchmark.sh
```

Na Windowsu koristite PowerShell ili Command Prompt:

```powershell
.\run-benchmark.ps1 -SuperLURoot C:\vcpkg\installed\x64-windows
```

```bat
run-benchmark.bat -SuperLURoot C:\vcpkg\installed\x64-windows
```

Skripte same ne preuzimaju zavisnosti; SuperLU i BLAS moraju biti već instalirani. Detaljna uputstva za macOS, Linux, Windows i WSL nalaze se u [docs/running.md](docs/running.md).

Podrazumijevani izlaz je jednostavna tabela sa četiri faze i ukupnim vremenom: vrijeme vaše implementacije, vrijeme SuperLU-a, razlika `ours - SuperLU (ms)` i odnos `ours / SuperLU (x)`. Pozitivna razlika znači da je vaša faza sporija; negativna znači da je vaša faza brža. Za dodatne dijagnostičke tabele koristite `--detailed` (Windows: `-Detailed`). CSV se upisuje samo ako eksplicitno proslijedite `--csv` (Windows: `-Csv`).

Ručno pokretanje na macOS-u i Linuxu:

```bash
make program
make bench
./benchmark/run.sh --repeat 3
```

Kompletan skup sa većim 2D matricama:

```bash
./benchmark/run.sh --repeat 3 --large --grid 50 --grid 60 --grid 80 --grid 100 --grid 120 --csv build/benchmark_results.csv
```

Benchmark mjeri svaku fazu odvojeno i podrazumijevano ispisuje samo vremena za četiri faze i ukupno vrijeme. Detaljna metodologija, svi rezultati i analiza uskog grla nalaze se u [docs/benchmark.md](docs/benchmark.md).

Glavni zaključak mjerenja na Apple M1: numerička Cholesky faktorizacija je konkurentna SuperLU-u (na istom poretku je čak brža u zidnom vremenu), dok je implementacija Minimum Degree orderinga dominantno usko grlo. Na 2D Laplacianu sa `n = 14400`, ordering zauzima oko 94% vremena i skoro 100% zaostajanja za SuperLU-om; popunjenost koju MD proizvodi je blizu SuperLU MMD popunjenosti, pa je problem u brzini implementacije MD-a (`std::set` i linearno traženje minimuma), a ne u kvalitetu poretka.

## Dokumentacija

Detaljnije objašnjenje cijelog procesa nalazi se u [docs/EXPLANATION.md](docs/EXPLANATION.md). Provjera ulazne matrice opisana je gore u sekciji "Provjera ulazne matrice".
