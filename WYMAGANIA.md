# Weryfikacja wymagań prowadzącego — OpenMP + CUDA

---

## 1. `default(none)` w `#pragma omp`

**Wymaganie:** Obowiązkowe stosowanie `default(none)` w pragmie OpenMP z jawnym wyspecyfikowaniem zmiennych prywatnych i współdzielonych.  
**Status:** ✅ OK (poprawione)

**Miejsce:** `OPENMP/openmp.c` (pętla główna)
```c
#pragma omp parallel for default(none) reduction(+:batch_found) schedule(dynamic) shared(bufory, graphs, hits, nhits, rf)
```

---

## 2. Rozdzielenie `#pragma omp parallel` i `#pragma omp for` przy alokacji pamięci

**Wymaganie:** Jeżeli w pętli alokowana jest pamięć lokalna wątków, rozbić pragmę na `#pragma omp parallel` + `#pragma omp for`.  
**Status:** ✅ OK

W pętli równoległej nie alokuje się pamięć — każdy wątek pracuje na już zaalokowanym buforze. Rozbicie niepotrzebne.

---

## 3. Unikanie `printf` w programach profilowanych

**Wymaganie:** Nie używać `printf` w pętli głównej podczas profilowania wydajności.  
**Status:** ✅ OK (poprawione)

Domyślnie oba programy działają w trybie cichym (quiet). Flaga `-d` włącza pełne raportowanie postępu:
```
./openmp n k [-a] [-d] [-b batch] [-m classes] [-B threads]
./cuda    n k [-a] [-d] [-b batch] [-m classes] [-B blockSize]
```

---

## 4. Wspólne interfejsy wejścia/wyjścia

**Wymaganie:** Dla obu programów rozwiązujących ten sam problem uzgodnić wspólne interfejsy.  
**Status:** ✅ OK

Oba programy mają identyczny interfejs wywołania:
```
./openmp n k [-a] [-d] [-b batch] [-m classes] [-B threads]
./cuda    n k [-a] [-d] [-b batch] [-m classes] [-B blockSize]
```
`-B` przyjmowane przez oba (OpenMP je ignoruje — brak bloków GPU). Wyniki zapisywane są do `result_n_k.txt` w tym samym formacie (graph6).

---

## 5. Unikanie konstrukcji angielskich w sprawozdaniu

**Wymaganie:** Nie używać anglicyzmów (np. "frameworku", "hotspot") w dokumentacji.  
**Status:** N/D (dotyczy sprawozdania, nie kodu)

---

## 6. Kompilacja w trybie Release (nie Debug)

**Wymaganie:** Do pomiarów czasowych kompilować z optymalizacjami `-O3`.  
**Status:** ✅ OK

Oba Makefile używają flag:
- OpenMP: `gcc -O3 -fopenmp`
- CUDA: `nvcc -O3 -arch=sm_86`

---

## 7. Unikanie `std::vector`

**Wymaganie:** Nie używać `std::vector` — utrudnia przeniesienie danych na GPU.  
**Status:** ✅ OK

Kod jest w czystym C — nie ma `std::vector`. Buforowane w surowych tablicach `malloc`.

---

## 8. Nie alokować/zwalniać wielokrotnie bez potrzeby

**Wymaganie:** Alokować raz na początku, zwalniać raz na końcu.  
**Status:** ✅ OK

- OpenMP: `malloc(batch_size * sizeof(char*))` + `malloc(BUFSIZE)` w pętli inicjalizacyjnej — jednorazowo przed pętlą główną. `free` na końcu.
- CUDA: `cudaMalloc` jednorazowo przed pętlą główną. `cudaFree` na końcu.

---

## 9. Niepotrzebne operacje w pętli `for`

**Wymaganie:** Przenosić stałe obliczenia przed pętlę.  
**Status:** ✅ OK

Główna pętla nie zawiera zbędnych powtórzeń. Wewnętrzny algorytm `eigensymmatrix` ma złożoną matematykę, ale operacje wewnątrz pętli są inherentną częścią algorytmu (nie można ich wynieść).

---

## 10. Parametry domyślne + możliwość zmiany z linii poleceń

**Wymaganie:** Nie zaszywać nazw plików w kodzie. Dane wejściowe najlepiej w formacie JSON. `int main()` zamiast `int main(int argc, char **argv)` to błąd.  
**Status:** ⚠️ CZĘŚCIOWO

- ✅ `int main(int argc, char **argv)` — OK.
- ✅ Parametry `n`, `k`, `-a`, `-b`, `-B` — OK.
- ❌ Nazwy plików `config_n_k.txt` i `result_n_k.txt` są zaszyte w kodzie. Brak flag `--config` / `--output`.
- ❌ Brak wejścia w formacie JSON.

---

## 11. Problemy grafowe — wiele instancji jednocześnie + kodowanie binarne

**Wymaganie:** Wersja CUDA powinna przetwarzać wiele instancji grafu jednocześnie. Struktura grafu kodowana binarnie. Mówimy o **liczbie** wierzchołków.  
**Status:** ✅ OK

- Każdy wątek GPU przetwarza jeden graf — batch_size instancji jednocześnie.
- Kodowanie graph6 jest binarne (6 bitów na bajt, spakowane w ASCII z offsetem 63).
- Terminologia: "liczba wierzchołków" (nie "ilość").

---

## 12. Nie stosować Pythona

**Wymaganie:** Python tylko za zgodą prowadzącego lub jako nakładka na bibliotekę C.  
**Status:** ✅ OK

Brak Pythona w projekcie.

---

## 13. CUDA — nie ustawiać sztywno `blockSize = 256`

**Wymaganie:** Pobrać `maxThreadsPerBlock` z GPU. Opcjonalnie dodać analizę czasową dla różnych rozmiarów bloków.  
**Status:** ✅ OK (poprawione)

`CUDA/cuda.cu` — `cudaGetDeviceProperties` pobiera `maxThreadsPerBlock` (na RTX 3060 = 1024). Flaga `-B` pozwala ręcznie nadpisać (z walidacją). Wcześniej było `int threadsPerBlock = 256` na sztywno — usunięte.

---

## 14. Funkcje biblioteczne thread-safe

**Wymaganie:** Nie używać `rand()` i innych funkcji niebezpiecznych wątkowo.  
**Status:** ✅ OK

- OpenMP: zapis do pliku chroniony `#pragma omp critical`.
- Obsługa sygnałów przez `sigaction` + `sigprocmask` jest thread-safe.
- Nie używamy `rand()`, `strtok()` itp.

---

## 15. Doxygen

**Wymaganie:** Dodać komentarze w konwencji Doxygen.  
**Status:** ❌ NIESPEŁNIONE

W kodzie nie ma ani jednego komentarza Doxygen (`@brief`, `@param`, `@return` itd.). Należy dodać dokumentację wszystkich funkcji.

---

## 16. Notacja matematyczna w sprawozdaniu

**Wymaganie:** $\mathcal{O}(n)$ zamiast $O(N)$. Wielkie litery oznaczają zbiory. Notacja jak w LaTeX.  
**Status:** N/D (dotyczy sprawozdania)

---

## 17. Składnia DokuWiki

**Wymaganie:** Zachować składnię DokuWiki.  
**Status:** N/D (dotyczy sprawozdania)

---

## 18. Generator liczb pseudolosowych

**Wymaganie:** Ziarno ustawiać w `main`, przekazywać jako parametr. Nie używać `rand()`.  
**Status:** N/D

Projekt nie używa generatorów pseudolosowych (grafy generowane przez `geng`).

---

## 19. Wznawianie wyszukiwania (geng res/mod)

**Wymaganie (metodyka zadania):** przestrzeń grafów dzielona klasami `geng ... res/mod`, praca wznawialna bez przeliczania od zera.  
**Status:** ✅ OK

- Oba programy przeszukują przestrzeń klasami: `geng -c n k:k r/mod` (domyślnie 256 klas, `-m` zmienia).
- `config_n_k.txt` przechowuje numery ukończonych klas (jedna liczba na linię, tryb append).
- Trafienia buforowane w pamięci osobno dla każdej klasy; zapis do `result_n_k.txt` tylko po pełnym przetworzeniu klasy.
- Przerwanie (SIGINT/SIGTERM) = odrzucenie bieżącej klasy; wznowienie kontynuuje od pierwszej nieskończonej klasy, bez duplikatów i bez pomijania po liczniku linii.

---

## Podsumowanie

| Status | Liczba |
|---|---|
| ✅ OK | 12 |
| ⚠️ Częściowo | 1 |
| ❌ Niespełnione | 1 |
| N/D | 5 |

### Do poprawy (w kodzie):

1. ~~**OpenMP `default(none)`**~~ — poprawione.
2. **Parametryzowane nazwy plików** — dodać flagi `--config` / `--output`.
3. ~~**Flaga `--quiet`**~~ — zrobione: domyślnie quiet, flaga `-d` włącza pełny output.
4. **Doxygen** — dodać komentarze dokumentacyjne.
