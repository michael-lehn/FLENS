# FLENS Build Compatibility and Regression Report

| | |
|---|---|
| Repository | https://github.com/michael-lehn/FLENS (branch `public`) |
| Baseline | commit `84043d90`, investigated 2026-10-26/27 (Part B) |
| Fixes | commits `f90daf7b` … `56731834`, validated 2026-10-27 (Part A) |
| Platform | macOS 26.6 (25G72), Apple M5, arm64. **No other platform was tested.** |
| Archive | harness, patches, reproductions, logs and raw results: `/Users/lehn/git/FLENS-compat/` (outside the repository) |

The report has two parts:

* **Part A** describes the current state: what was fixed (with commits), the final validation, corrections to the original findings, issues found while fixing, and the remaining limitations.
* **Part B** is the original investigation of the unmodified commit `84043d90`. It is kept as it was, apart from the corrections marked **[Correction]** and the status column in §8. Section numbers in Part A (§6.x, §7.x) refer to Part B.

---

## 0. Executive summary

**Original state (`84043d90`):**
* FLENS worked with modern compilers only in C++11/C++14 mode. C++20 and later failed with every compiler, and C++17 failed with libc++.
* Several genuine numerical bugs existed independently of the compiler. The most serious: `lapack::ev(false, HeMatrix, …)` / `zheev('N')` returned wrong eigenvalues.
* Two Fortran-callable LAPACK routines crashed (`zheev_` with N=0, `zgelqf_` with small workspace).
* The build tooling did not work out of the box: hardcoded `g++-4.8`, mandatory `-Werror`, missing `check` targets, no header dependency tracking, CMake 4 rejected, and gfortran ≥ 10 rejected the Netlib sources.

**Current state (`56731834`, tested on macOS/arm64 only):**
* **All five numerical bugs that were reachable** (sterf loop, sterf `SafeMin`, tiny `copy`/`ccopy`, tiny `gemv`, `lan`) **and both LAPACK crashes are fixed.** Each fix has a regression test that fails before the fix and passes after it, and is checked against reference LAPACK where one exists. The sixth reported bug (`jsv.tcc` `Eps`) turned out to be in commented-out code; see A.3.
* Two further bugs found during the work are fixed: `lsame_` in `lapack_FLENS.a` compared case-sensitively, and the MatrixMarket loaders aborted under most locales with libstdc++.
* **C++11, 14, 17, 20, 23 and 26 now give identical results** on GCC 14/15/16, Clang 17/20/23 and Apple Clang 21: the same set of programs passes, and stdout is bit-identical across standards. C++11 still works.
* **Netlib LAPACK 3.3.1, extended d/z set:** FLENS-LAPACK now gives **exactly the same results as pure reference LAPACK** in every test that FLENS supports (1,831,162 tests). The only exception is dcsd/zcsd (540 tests), which needs a workspace query that FLENS does not implement.
* **`dev-check` passes bitwise:** 956,515 tests, each FLENS-LAPACK result compared bit for bit with reference LAPACK.
* **`make check`, `check-blas`, `check-lapack`, `check-lapack-extended` and `dev-check` work** with configurable compilers, evaluate the results and fail on unexpected results. `-Werror` is opt-in, header changes trigger rebuilds, and CMake 4 is supported.
* **Not fixed (out of scope or a design decision):** the C++-standard-independent compile failures of some examples, tutorials and playground programs (§6.2.10–6.2.17, §6.4), the `dgejsv_` interface, the `orgqr`/`ungqr` workspace query, and the biased-pointer UB pattern. Details in A.5.

---

# Part A — Fixes and final validation

## A.1 Commits

All commits are on `public` and pushed. Each fix has its own commit, and every fix commit adds a regression test in `flens/test/regression/`.

| Commit | Change | Finding | Regression test |
|---|---|---|---|
| `f90daf7b` | `sterf`: QL search loop (`m` shadowed, `l` incremented, `m` unset for `l==lEnd`) | §6.2.3 | `sterf.cc` |
| `2131b158` | `sterf`: `One / SafeMin` → `One / safeMin` | §6.2.4 | `sterf.cc` (wide-range cases) |
| `3b9d5be7` | tiny `copy` and `ccopy`: `y[iX]` → `y[iY]` | §6.2.1 (+ `ccopy`, same bug) | `tiny-copy.cc` |
| `5dd06ec5` | tiny `gemv`: do not read y for beta = 0; `iX += incX` in `gemv_t`/`gemv_ct` | §6.2.2 | `tiny-gemv.cc` |
| `2da9dc2a` | `lan`: `++j` → `++i` (infinity norm, upper non-unit `TrMatrix`) | §6.2.8 | `lantr.cc` |
| `6d224bd8` | `zheev_`: RWORK view length `max(1, 3N-2)` | §6.2.6 | `zheev-interface.cc` |
| `98a72f89` | `lqf`: `""` instead of a null pointer for the `ilaenv` option string | §6.2.7 | `lqf-ilaenv.cc` |
| `909570db` | MatrixMarket loaders: `std::locale::classic()` instead of `std::locale("")` | new (A.4) | `matrixmarket.cc` |
| `9fffd7f6` | C++17/C++20: `std::allocator_traits`; lambda instead of `bind2nd`/`ptr_fun` | §6.1.1, §6.1.2 | whole matrix (A.2) |
| `14878d34` | `lsame_` in `lapack_FLENS.a`: case-insensitive | new (A.4) | `lsame-interface.cc` |
| `17daf8e5` | `toolchain.mk`: configurable `CXX`/`CXXSTD`/`CXXFLAGS`/`FC`, `WERROR=1` opt-in, `-MMD -MP`, fixed `make dev` archive | §6.3.1, 6.3.3, 6.3.4 | build checks (A.2.5) |
| `7914b0a1` | Netlib reference BLAS/LAPACK: `-fallow-argument-mismatch` (probed) and `-ffp-contract=off` | §6.3.7 | Netlib suites |
| `5a710c11` | `make check`, `check-regression`, `check-blas`, `check-lapack`, `check-lapack-extended`, `dev-check`; results evaluated (`check-results.sh`, `check-expected.txt`) | §6.3.2 | — |
| `9529905e` | CMake: `cmake_minimum_required(VERSION 3.5...3.31)` before `project()` | §6.3.6 | CMake 4.4 configure + install |
| `56731834` | `.gitignore`: unpacked tarballs, dependency files | — | — |

Regression tests: `make check-regression`, or `make -C flens/test/regression check`. `zheev-interface` and `lqf-ilaenv` need a LAPACK library (`REFLAPACK=...`); with it, the other tests are also compared with reference LAPACK. Each test was confirmed to **fail on the parent commit and pass with the fix**.

**Searches for similar mistakes** (done for each fix):
* Enumerators used as values in arithmetic, such as `SafeMin`/`Eps`: one more hit, `jsv.tcc:907`, but it is dead code (A.3).
* `for` loops whose increment does not match the loop variable (pattern search over `flens`, `cxxblas`, `cxxlapack`; it finds the `lan.tcc` bug in the baseline): no other hits.
* Destination indexed with the source index in the tiny kernels: `ccopy` had the same bug and was fixed in the same commit.
* Kernels that read y when beta = 0: only the tiny `gemv`. The generic dense `gemv`/`gemm` use `scal(beta, y)` (§8 item 24, unchanged).
* `ilaenv` called with `0` as the option string: no other hits.
* Interface wrappers creating views of length `3N-2` or similar: no other hits.
* `std::locale("")`: only the two MatrixMarket loaders.

## A.2 Final validation (HEAD `56731834`)

### A.2.1 Compiler × standard matrix (`-O2`)

Same harness and program inventory as Part B §4: 7 compilers × 6 standards × 129 programs, 5,418 builds. `test/geccsmv` is excluded and run separately (A.2.2).

| Compiler | C++11 | C++14 | C++17 | C++20 | C++23 | C++26 |
|---|---|---|---|---|---|---|
| gcc14, gcc15, gcc16 | 96 / 27 / 2 | same | same | same | same | same |
| clang17 | 98 / 27 / 0 | same | same | same | same | same |
| clang20, clang23, appleclang21 | 97 / 28 / 0 | same | same | same | same | same |

Columns are ran & passed / compile error / link error. In every configuration: 0 timeouts, 0 crashes, 3 expected compile failures confirmed, 1 skipped (MKL).

**Compared with the original baseline (§4.1):**
* The **C++11 counts are unchanged.** The baseline's only timeout was `geccsmv`, which is now run separately.
* **C++17, C++20, C++23 and C++26 went from 1–4 passing programs to the same result as C++11** on all 7 compilers.
* The Clang C++26 failure of `lapack-complex-heev` seen with the minimal patch alone (§4.2) is gone, because the `SafeMin` fix removed the enum/float arithmetic.
* The remaining 31–33 non-passing programs per compiler are exactly those already analysed in §6.2.10–6.2.15 and §6.4. They do not depend on the C++ standard. The difference between compilers comes from Clang 17's handling of the MPFR hack (§6.2.13) and from §6.2.11, which affects libstdc++ only.

**`-O0`** (GCC 16 and Clang 23, C++20): same counts as at `-O2`.

**ASan + UBSan** (GCC 16 and Clang 23, C++20, `-O1 -g`): every program has the same status as in the C++11 baseline run (§7.2): 89 and 69 clean runs. The UBSan reports are the same biased-pointer pattern on empty views. There are no ASan reports.

### A.2.2 Long-running test

`test/geccsmv` (3000 random runs), GCC 16, C++20, `-O2`: **passed** (1,143 s).

### A.2.3 Numerical output compared with the baseline

The stdout of every program was compared with the original baseline outputs (`compare_outputs.py`):
* **GCC 16:** 96 of 96 bit-identical at C++11, C++20 and C++26.
* **Clang 17/23, Apple Clang 21:** all bit-identical except **`examples/tiny_gemv`**, which now prints the correct `10 26 30 32` instead of `nan` (fix `5dd06ec5`). Same at C++11, C++20 and C++26.
* Across C++ standards, outputs are bit-identical for each compiler.
* No other example output changed. None of the shipped examples calls `ev(false, HeMatrix)`, `sterf` on wide-range data, tiny copies with different strides, or `lan` on upper triangular matrices. The effects of those fixes show up only in the regression tests and in the Netlib suites.

### A.2.4 Netlib test suites

**BLAS** (`make check-blas`): 40 of 40 level 1, 96 of 96 level 2/3, 78 of 78 error exits. Built with the default toolchain, GCC 16 C++11/C++20, and Clang 23 C++20.

**LAPACK 3.3.1, `make check-lapack`** (zsg, dtest, ztest, dnep, znep): 956,515 tests, 0 failures. Built with the default toolchain (Apple Clang, C++11), GCC 16 C++11/C++20, and Clang 23 C++20.

**`make dev-check`** (FLENS-LAPACK compares each result bit for bit with reference LAPACK): 956,515 tests, 0 failures with GCC 16 C++11 and Clang 23 C++20. This needs `-ffp-contract=off` for both the C++ and the Fortran side, plus `-fcx-fortran-rules` for complex division. `toolchain.mk` adds both automatically when the compilers accept them (`DEV_CXXFLAGS`, `FC_REF`). Without them, dev-check reports differences in the last bits.

**Extended d/z EIG set** (`make check-lapack-extended`; GCC 16 C++11 via make; GCC 16 C++20 and Clang 23 C++20 via the harness):

| | Baseline `84043d90` | After fixes | Pure reference LAPACK 3.3.1 |
|---|---|---|---|
| Tests run (completed blocks) | 1,786,520 | 1,831,162 | 1,831,702 |
| zsep | **abort** (`zheev_`, N=0) | 3 threshold misses (same as reference) | 3 threshold misses |
| zsvd | **segfault** (`zgelqf_`) | 1 threshold miss (same as reference) | 1 threshold miss |
| zgbak | passed, but ZGEMM rejected `'n'` as an illegal value (`lsame_`, A.4) | clean | clean |
| dsep, dgd, zgd | threshold misses | same as reference | same |
| dcsd, zcsd | abort (`orgqr`/`ungqr` workspace query) | **abort (unchanged)** | 270 + 270 tests pass |
| all other d/z inputs | pass | pass | pass |

* The 540-test difference from reference LAPACK is exactly dcsd + zcsd.
* `check-results.sh` reports **32 passed, 8 expected deviations, 0 failed**. The 8 expected deviations are listed exactly in `flens/lapack/check-expected.txt`:
  * 5 are threshold misses identical to reference LAPACK (dgd, zgd, dsep, zsep, zsvd).
  * 2 are dcsd/zcsd.
  * 1 is dsvd's DGEJSV message (A.3).
* The per-test summaries (test counts, threshold misses) after the numerical fixes (`98a72f89`), after the compatibility changes (`14878d34`, GCC 16 and Clang 23 at C++20) and at HEAD (GCC 16, C++11) are identical.

### A.2.5 Build system

| Check | Result |
|---|---|
| `make check` with defaults (`c++` = Apple Clang 21, `gfortran`) | passes (regression, BLAS, LAPACK) |
| `make check CXX=g++-16`, `CXX=clang++` (LLVM 23) `CXXSTD=c++20` | passes |
| `make -j check dev-check` | passes; the LAPACK targets are serialized because they share `flens/lapack/test` |
| No-op rebuild | nothing rebuilt |
| `touch` a header | only the dependent objects are rebuilt (`-MMD -MP`) |
| `make distclean` | the tree is pristine again (compared with `git status --ignored`) |
| `make -j` (libraries, Netlib reference builds) | no races observed |
| `WERROR=1` | **fails**: real warnings remain (A.5) |
| CMake 4.4: configure, build, `cmake --install` | works; a tutorial compiles against the installed headers |

### A.2.6 Regression tests across compilers and standards

All 7 compilers × C++11/14/17/20/23/26 (42 configurations), with `REFLAPACK` set to OpenBLAS: **8 of 8 tests pass in each configuration**.

## A.3 Corrections to the original report

* **§6.2.5 / §8 item 6 (`jsv.tcc:907`, `Eps`): the line is inside commented-out code.** It is not compiled, so it is not a live bug and was not changed. Investigating it showed that `jsv` as a whole does not work:
  * The generic `lapack::jsv` is a stub (`ASSERT(0)`).
  * The interface `dgejsv_` (`flens/lapack/interface/src/gejsv.cc:35–40`) has a broken LWORK check (`||` where `&&` is needed). It therefore always returns INFO = −17, so `dgejsv_` never computes anything.
  * The Netlib driver `ddrvbd` checks the stale `IINFO` instead of `INFO`, which is why dsvd still counts as passing. The message "On entry to DGEJSV parameter number 17 had an illegal value" is listed as an expected deviation.
  * Making `jsv` work means porting `DGEJSV`, not fixing one line. That is a design decision, so it was left open.
* **§6.2.4 (`SafeMin`): the effect was understated.** The original report said the scaling is undone afterwards. In fact, blocks with a large norm are scaled down until their small entries become subnormal or zero, which gives wrong eigenvalues. For example, diag(1e200, 1, 1) with weak coupling gave {1, 1, 1e200} instead of {0, 2, 1e200}. This is now covered by the regression test.
* **§5.3: zgbak was counted as passing, but `lsame_` was broken.** The zgbak output contains binary bytes, and `grep` treated it as a binary file, so the ZGEMM "illegal value" messages were missed. The evaluator now uses `grep -a`.
* **§6.1.2 suggested `std::locale("")` in the lambda.** The committed version uses `std::locale::classic()` (A.4).

## A.4 Issues found during the work

**Fixed:**
* **`lsame_` in `lapack_FLENS.a` was case-sensitive** (`14878d34`). Reference LAPACK routines linked with `lapack_FLENS.a` therefore rejected lower-case options.
* **MatrixMarket loaders used `std::locale("")`** (`909570db`). With libstdc++ on macOS this throws for every locale except C/POSIX, so loading a file aborted the program, for example under `LANG=en_US.UTF-8`.

**Not fixed (documented here):**
* MatrixMarket loading into `SyMatrix`/`HeMatrix` calls `resize(numRows, numCols)`. For these types the arguments mean `(dim, firstIndex)`, so loading hits an index assertion.
* Packed MatrixMarket save/load for `SpMatrix`/`HpMatrix`, and load for `HpMatrix`, do not compile when instantiated (they use `firstIndex`/`lastIndex`).
* `gecopy` with `ConjTrans` does not conjugate; there is also a `trans==Conj || trans==Conj` typo.
* `lasrt`: its biased stack pointer triggers GCC `-Warray-bounds`. A `-Wmaybe-uninitialized` warning for `sigma` is a false positive.
* Clang 23 warns `-Wuninitialized-const-pointer` for `DUMMY` in `flens/lapack/impl/ormqr.tcc:291` and `ormlq.tcc:289`. These are workspace-query dummies that are not read.

## A.5 Remaining limitations

| Item | Classification | Status |
|---|---|---|
| dcsd/zcsd abort: `orgqr`/`ungqr` workspace query (`LWORK=-1`) not implemented (§6.2.9) | FLENS limitation | expected deviation in `check-lapack-extended` |
| `dgejsv_` always returns INFO = −17; generic `jsv` is a stub (A.3) | FLENS limitation / design decision | expected deviation (dsvd) |
| Threshold misses in dgd, zgd, dsep, zsep, zsvd | test suite / toolchain; identical with reference LAPACK | expected deviation |
| `WERROR=1` fails: `-Wunused-local-typedefs`, `-Warray-bounds` (`lasrt`), `-Wmaybe-uninitialized`, Clang `-Wuninitialized-const-pointer` | pre-existing warnings | `-Werror` is opt-in; default builds are clean of errors |
| `FLENS_DEBUG_CLOSURES` broken (`x1x2`, §6.2.10): 3 tutorials | FLENS bug, compile error, all standards | not fixed (outside the requested scope) |
| `std::complex<int>` from index-variable closures, libstdc++ only (§6.2.11): 2 examples | FLENS bug, compile error | not fixed |
| `tut08-page01-example2/5` stale API (§6.2.12) | example bug | not fixed |
| MPFR/QD trait specializations (§6.2.13) | FLENS UB; compile errors on Clang ≥ 20 | not fixed |
| Playground: `cgs.h` include guard (§6.2.14), `ptfqmr`, `suitesparse`, SuperLU header nesting (§6.2.15) — 18 programs | playground bugs | not fixed |
| `gecrs-superlu` (SuperLU ≥ 5 API), QD + libstdc++, MPI C++ bindings, MKL (§6.4) | external dependencies / environment | not fixed |
| Biased-pointer UB on empty views (§7.2) | latent UB, deliberate design | not fixed; no wrong results observed |
| Generic `gemv`/`gemm` propagate NaN from y for beta = 0 | deviation from BLAS semantics | not fixed; not observable with zero-initialized dense storage |
| Issues in A.4 ("not fixed") and §6.2.16, §6.2.17 | latent | not fixed |

**Compatibility, tested configurations only:**
* macOS 26.6 on arm64.
* GCC 14.4/15.3/16.2 with libstdc++, Clang 17/20/23 with Homebrew libc++, Apple Clang 21 with the SDK libc++.
* gfortran 16.2; GNU Make 4.4.1; CMake 4.4.4.
* Standards: C++11 to C++26.

**Not tested** (see also §10): Linux, x86_64, Windows/MSVC, Intel compilers, older GCC/Clang releases, other Fortran compilers, MKL, and running the MPI programs. No compatibility claim is made for these.

---

# Part B — Original investigation of `84043d90`

The following sections describe the unmodified baseline. Paths in §9 now point to the permanent archive.

## 1. Environment and compiler versions

| Tool | Version | Path |
|---|---|---|
| GCC | 14.4.0 (Homebrew `gcc@14`) | `/opt/homebrew/bin/g++-14` |
| GCC | 15.3.0 (Homebrew `gcc@15`) | `/opt/homebrew/bin/g++-15` |
| GCC | 16.2.0 (Homebrew `gcc`) | `/opt/homebrew/bin/g++-16` |
| Clang | 17.0.6 (Homebrew `llvm@17`) | `/opt/homebrew/opt/llvm@17/bin/clang++` |
| Clang | 20.1.8 (Homebrew `llvm@20`) | `/opt/homebrew/opt/llvm@20/bin/clang++` |
| Clang | 23.1.4 (Homebrew `llvm`) | `/opt/homebrew/opt/llvm/bin/clang++` |
| Apple Clang | 21.0.0 (clang-2100.3.34.2) | `/usr/bin/clang++` |
| gfortran | 16.2.0 (Homebrew GCC) | `gfortran-16` |
| GNU Make | 4.4.1 | |
| CMake | 4.4.4 | |
| Python | 3.15.0 (harness only) | |

* **Standard libraries:**
  * The GCC builds use libstdc++.
  * The Homebrew Clang builds use Homebrew libc++.
  * Apple Clang uses the SDK's libc++.
* **C++26:** every compiler accepted `-std=c++26`, so no `c++2c` fallback was needed.
* **System compiler configuration:** not changed. All compilers were invoked by absolute path.

---

## 2. Dependencies and installation steps

Already present: Homebrew GCC 16, LLVM 17–23, Boost 1.92.0, GMP 6.3.0, MPFR 4.2.2, CMake, GNU Make.

Installed for this investigation:

```sh
HOMEBREW_NO_AUTO_UPDATE=1 brew install openblas lapack superlu qd fftw open-mpi suite-sparse
HOMEBREW_NO_AUTO_UPDATE=1 brew install gcc@15 gcc@14
# MTL4 (header-only; needed only by examples/design-flaws-mtl4.cc), scratch dir only:
git clone --depth 1 https://github.com/simunova/mtl4.git <scratch>/deps/mtl4   # commit 09d8523 (2021-03-05)
```

| Package | Version | Used for |
|---|---|---|
| openblas | 0.3.34 | LAPACK backend for `-DUSE_CXXLAPACK` examples (replaces the docs' `-framework vecLib`); CBLAS backend experiment |
| Apple Accelerate | macOS 26.6 SDK | docs-faithful vecLib path (Apple Clang, LAPACK examples) |
| lapack | 3.12.1 | installed but not needed in the end (reference LAPACK 3.3.1 comes from the repo's own `cxxlapack/netlib/lapack.tgz`) |
| superlu | 7.0.1 | `geccs-superlu`, `gecrs-superlu`, `playground/example-superlu` |
| qd | 2.3.24 | `qd-lapack-*` |
| mpfr / gmp | 4.2.2 / 6.3.0 | `mpfr-real-lapack-*` |
| fftw | 3.3.11 | `playground/example-fft-*` |
| open-mpi | 5.0.11 | `playground/example-mpi-*` |
| suite-sparse | 7.14.1 | `playground/example-suitesparse` |
| boost | 1.92.0 | MTL4 prerequisite |
| MTL4 | git 09d8523 | `design-flaws-mtl4` |

Intel MKL (`playground/example-mkldss`) is not available for macOS/arm64 and was **skipped**.

---

## 3. Inventory of examples and tests

### 3.1 Build systems found

| Item | What it does | Status |
|---|---|---|
| top-level `Makefile` | `all` prints "header only" hint; `check-blas`, `check-lapack`, `dev-check`, `clean`, `distclean` | `make check` (named in README) **does not exist**; `check-lapack` → `make -C flens/lapack check` → *No rule to make target 'check'*; `dev-check` → `make -C flens/blas/interface dev-check` → *No rule*. |
| `CMakeLists.txt` (+ per-dir) | install headers only | fails with CMake ≥ 4 (see §6.3) |
| `cxxblas/netlib/Makefile` | unpacks + builds Netlib reference BLAS/CBLAS | works (also with `-j`) |
| `cxxlapack/netlib/Makefile` | unpacks + builds reference LAPACK 3.3.1 | needs `-fallow-argument-mismatch` with gfortran ≥ 10 |
| `flens/blas/interface/{,blas/}Makefile` | builds `blas_FLENS.a` (Fortran-callable FLENS-BLAS) + Netlib BLAS tests | hardcodes `CXX = g++-4.8`, `-Werror` |
| `flens/lapack/{,interface/}Makefile` | builds `lapack_FLENS.a` + Netlib LAPACK 3.3.1 test suite (`dev-check`) | hardcodes `CXX = g++-4.8`, `-Werror`; `dev` target archives the wrong objects |

**No Makefile builds or runs the examples, `flens/test`, or `playground`.** The only build instructions for these are the compile commands written in the `*.doc` tutorial files. They are therefore all "not covered by the default build". Every one of them was compiled and run by the harness described in §9.

**Required standard:** the tutorials compile with `-std=c++11`, and the README says `g++ (>= 4.7)`. Nothing claims support for C++14 or later. All results above C++11 are **experimental observations, not intended support**.

### 3.2 Programs with `main()` — 130 total

| Group | Count | Notes |
|---|---|---|
| `flens/examples/*.cc` | 103 | tutorials (`tut01`…`tut08`), LAPACK examples (`lapack-*`, 56), sparse, tiny, MPFR/QD, SuperLU, MTL4 comparison |
| `flens/test/*.cc` | 8 | self-checking sparse mat-vec tests (`ASSERT`-based; compare sparse vs. dense results; 30 random runs; `geccsmv` runs **3000**) |
| `playground/flens/examples/*.cc` | 19 | Krylov solvers, FFT (FFTW), MPI, SuperLU, SuiteSparse, MKL |

**Build profiles taken from the `.doc` files:**
* `lapack-*` and `tut05-*`: `-DUSE_CXXLAPACK`, linked against OpenBLAS LAPACK.
* `tut02-page03/04`: `-DFLENS_DEBUG_CLOSURES` plus `flens/debug/auxiliary/*.cc`.
* SuperLU, QD and MPFR examples: their respective libraries.
* Playground programs: libraries chosen by the `#define WITH_*` in each source.
* Everything else: header-only (`-I<root>`).

**Intended compile failures** (the tutorial text says they must not compile):
* `tut01-page03-example04`
* `tut05-page02-example2`
* `tut05-page02-example4`

These count as passed when they fail with the documented diagnostic.

**Special cases (documented intent, verified):**
* `vecadd` returns `y(2)` as its exit status. Exit code **19** = 2·1 + 4·2 + 3·3 is the correct result.
* `addingvectors` and `addingvectors_bench` abort by design without `-DUSE_SIMPLE_ET`. FLENS asserts on `y += x1 + x2`, which would need a temporary. The commit introducing them ("alternate evaluation strategies") points to `simple_et.h`, so they are run with `-DUSE_SIMPLE_ET`. The plain-build aborts were also recorded (SIGABRT on all 7 compilers).
* `playground/example-fft-import-wisdom` consumes the `wisdom.dat` written by `example-fft-export-wisdom`.
* MPI programs are run with `mpirun -np 2 --oversubscribe`.

### 3.3 Validation suites

* **Netlib BLAS test suite** (`flens/blas/interface/test.tgz`):
  * `[sdcz]blat1`: 40 subprograms.
  * `[sdcz]blat2/3`: 96 computational + 78 error-exit checks.
* **Netlib LAPACK 3.3.1 test suite** (`flens/lapack/test.tgz`):
  * `dev-check` runs zsg, dtest, ztest, dnep, znep.
  * Additionally run here: **all 38 d/z EIG test inputs**.
  * Only d/z were tested, because `lapack_FLENS.a` exports only double and complex-double routines (31 `d*` and 25 `z*` symbols).

---

## 4. Build compatibility matrix

The tables count **all 130 programs** for each compiler and standard (`-O2`, unmodified sources). "ran & passed" means: compiled, linked, exited with the expected status, printed no assertion text, and finished within the timeout. Numerical validation is covered in §5 and §7.

### 4.1 Unmodified repository (`-O2`)

| Compiler | Std | ran & passed | compile error | link error | timeout¹ | expected-fail confirmed | expected-fail, other reason² | skipped³ |
|---|---|---|---|---|---|---|---|---|
| gcc14 | C++11 | 96 | 27 | 2 | 1 | 3 | 0 | 1 |
| gcc14 | C++14 | 96 | 27 | 2 | 1 | 3 | 0 | 1 |
| gcc14 | C++17 | 96 | 27 | 2 | 1 | 3 | 0 | 1 |
| gcc14 | C++20 | 4 | 122 | 0 | 0 | 3 | 0 | 1 |
| gcc14 | C++23 | 4 | 122 | 0 | 0 | 3 | 0 | 1 |
| gcc14 | C++26 | 4 | 122 | 0 | 0 | 3 | 0 | 1 |
| gcc15 | C++11/14/17 | 96 | 27 | 2 | 1 | 3 | 0 | 1 |
| gcc15 | C++20/23/26 | 4 | 122 | 0 | 0 | 3 | 0 | 1 |
| gcc16 | C++11/14/17 | 96 | 27 | 2 | 1 | 3 | 0 | 1 |
| gcc16 | C++20/23/26 | 4 | 122 | 0 | 0 | 3 | 0 | 1 |
| clang17 | C++11/14 | 98 | 27 | 0 | 1 | 3 | 0 | 1 |
| clang17 | C++17/20/23/26 | 1 | 125 | 0 | 0 | 3 | 0 | 1 |
| clang20 | C++11/14 | 97 | 28 | 0 | 1 | 3 | 0 | 1 |
| clang20 | C++17/20 | 1 | 125 | 0 | 0 | 3 | 0 | 1 |
| clang20 | C++23/26 | 1 | 125 | 0 | 0 | 2 | 1 | 1 |
| clang23 | C++11/14 | 97⁴ | 28 | 0 | 1 | 3 | 0 | 1 |
| clang23 | C++17/20 | 1 | 125 | 0 | 0 | 3 | 0 | 1 |
| clang23 | C++23/26 | 1 | 125 | 0 | 0 | 2 | 1 | 1 |
| appleclang21 | C++11/14 | 97⁴ | 28 | 0 | 1 | 3 | 0 | 1 |
| appleclang21 | C++17/20 | 1 | 125 | 0 | 0 | 3 | 0 | 1 |
| appleclang21 | C++23/26 | 1 | 125 | 0 | 0 | 2 | 1 | 1 |

(Rows merged where the per-standard counts are identical. 42 configurations × 130 programs = 5,460 records.)

¹ **Timeout** is `test/geccsmv`, which runs 3000 random cases where its siblings run 30. It was re-run without a time limit and **passed all 3000 runs**: GCC 16 in 2364 s, Clang 23 in 2110 s.

² In C++17 and later every program already fails on the C++17/20 blockers. The "expected compile failures" therefore fail for the wrong reason, and their "confirmed" status above C++14 is meaningless.

³ `playground/example-mkldss` needs Intel MKL.

⁴ Includes `examples/tiny_gemv`, which exits 0 but prints **`nan`** under every Clang at `-O2` (UB, §6.2.2). Numerically valid count: 96.

The few "passes" above C++17 are programs that never instantiate the broken code:
* `design-flaws-mtl4` doesn't use FLENS at all.
* With GCC in C++20 and later, `tinyarray`, `tinyfullstorage` and `tiny_gemv` also pass, since tiny storage uses no allocator.

**C++11 status per group** (GCC 16 / Clang 23):
* examples: 92 / 93 of 103 pass (includes the 3 expected failures).
* test: 7 of 8 pass within the timeout, 8 of 8 when given unlimited time.
* playground: **0 of 19** (§6.3.4).

### 4.2 Same matrix with only the minimal C++17/20 patch applied (scratch copy)

The patch file is listed in §9. GCC 16, Clang 23 and Apple Clang 21 were tested across all 6 standards; `test/geccsmv` was excluded for time (2,304 builds).

| Compiler | C++11 | C++14 | C++17 | C++20 | C++23 | C++26 |
|---|---|---|---|---|---|---|
| gcc16 (ran & passed / compile err / link err) | 96/27/2 | 96/27/2 | 96/27/2 | 96/27/2 | 96/27/2 | 96/27/2 |
| clang23 | 97/28/0 | 97/28/0 | 97/28/0 | 97/28/0 | 97/28/0 | **96/29/0** |
| appleclang21 | 97/28/0 | 97/28/0 | 97/28/0 | 97/28/0 | 97/28/0 | **96/29/0** |

* The 3 expected failures are confirmed in every configuration.
* With the patch, **C++17 through C++23 are identical to C++11**: same pass and fail sets, and **bit-identical stdout** for all 96/97 programs.
* The only new failure is `lapack-complex-heev` under Clang in C++26. C++26 makes arithmetic between a floating-point value and an enumeration ill-formed (P2864), and that exposes the real `sterf` typo (§6.2.4). GCC 16 accepts it in C++26 with a `-Wdeprecated-enum-float-conversion` warning.

### 4.3 Summary of intended vs. observed support

| Standard | Intended (docs) | Observed, unmodified | Observed, with 16-line patch |
|---|---|---|---|
| C++11 | yes | works on all 7 compilers (except the C++11-independent defects in §6) | works |
| C++14 | not stated | identical to C++11 | identical to C++11 |
| C++17 | not stated | GCC/libstdc++: works (deprecation warnings). libc++ (all Clang): **every program fails** | identical to C++11 |
| C++20 | not stated | **every program using dense storage fails on every compiler** | identical to C++11 |
| C++23 | not stated | fails | identical to C++11 |
| C++26 | not stated | fails | identical to C++11, except one Clang error that exposes a real bug |

---

## 5. Execution results

### 5.1 Totals

| Campaign | Builds | Program executions |
|---|---|---|
| Baseline matrix, 7 compilers × 6 standards, `-O2` | 5,418 + 210 re-runs⁵ | 1,711 |
| Minimal-patch matrix, 3 compilers × 6 standards | 2,304 | 1,738 |
| `-O0`, 7 compilers, C++11 | 896 | 677 |
| ASan+UBSan, GCC 16 + Clang 23, C++11, `-O1 -g` | 256 | 193 |
| Strong warnings (compile only), GCC 16 + Clang 23 | 258 | – |
| Alternative backends (OpenBLAS CBLAS, reference LAPACK, Accelerate) | 222 | 204 |
| Playground with include-guard fix | 36 | 24 |
| **Total (harness)** | **9,390** | **4,547** |

⁵ The re-runs repeated 5 programs over all 42 configurations after their documented intent was established (§3.2).

Additionally: 2 unlimited-time `geccsmv` runs, plus the Netlib suites below.

### 5.2 Netlib BLAS test suite against FLENS-BLAS (`blas_FLENS.a`)

| Build of `blas_FLENS.a` | Level 1 (40) | Level 2/3 computational (96) | Error exits (78) | Failures |
|---|---|---|---|---|
| GCC 16, C++11, default flags (no `-O`) | 40 PASS | 96 PASSED | 78 PASSED | **0** |
| GCC 16, C++11, `-O2` | 40 | 96 | 78 | **0** |
| Clang 23, C++11, `-O2`, libc++ linked into gfortran | 40 | 96 | 78 | **0** |
| GCC 16, **C++20**, `-O2`, minimal-patch copy | 40 | 96 | 78 | **0** |

`CBLAT3`/`ZBLAT3` themselves print "ERROR-EXITS WILL NOT BE TESTED"; that is a setting in the shipped `.dat` files.

### 5.3 Netlib LAPACK 3.3.1 test suite against FLENS-LAPACK (`lapack_FLENS.a`, GCC 16, C++11)

**`dev-check` set (the project's own check target):**

| Test | Tests run | Failures |
|---|---|---|
| zsg | 30,870 | 0 |
| dtest (LIN) | 572,610 | 0 (+ all error-exit tests passed) |
| ztest (LIN) | 325,105 | 0 (+ all error-exit tests passed) |
| dnep | 13,230 | 0 |
| znep | 14,700 | 0 |
| **Total** | **956,515** | **0** |

**All other d/z EIG inputs** (not part of `dev-check`; run here to exercise more FLENS code). Every failure was re-run against **pure reference LAPACK** (FLENS archive replaced by an empty one) to attribute it.

| Test | FLENS-LAPACK | Pure reference LAPACK 3.3.1 | Attribution |
|---|---|---|---|
| dbb, dec, ded, dgg, dglm, dgqr, dgsv, dlse, dsb, dsg, dsvd, dbak, dbal, dgbak, dgbal | all pass (6,000 / 501,251 / 20,080 / 13,832 / 48 / 1,728 / 384 / 96 / 540 / 30,870 / 73,000 tests) | – | – |
| zbb, zec, zed, zgg, zglm, zgqr, zgsv, zlse, zsb, zbak, zbal, zgbak, zgbal | all pass¹ | – | – |
| dsep | 1/4662 DST + 1/14256 DST drivers over threshold | **same** | test suite / toolchain |
| dgd | 200/5000 DXV over threshold | **same** | test suite / toolchain |
| zgd | 24/5000 ZXV over threshold | **same** | test suite / toolchain |
| **zsep** | **abort**: `Assertion (length_>=0)` in `ArrayView` via `zheev_` | 1+1+1 over threshold, completes | **FLENS bug** (§6.2.5, §6.2.6) |
| **zsvd** | **segfault** in `strlen` via `zgelqf_` → `ilaenv` | 1/4300 ZBD over threshold, completes | **FLENS bug** (§6.2.7) |
| **dcsd / zcsd** | **abort**: `ASSERT(0)` in `dorgqr_`/`zungqr_` | 270 tests, all pass | **FLENS limitation**: workspace query is a `// TODO` |

¹ **[Correction]** zgbak printed "illegal value" messages from ZGEMM, caused by a case-sensitive `lsame_` in `lapack_FLENS.a`. The evaluator missed them because the output file contains binary bytes (A.3, A.4).

**Verified in a scratch copy** (patch in §9):
* With the `heev` (N=0) fix alone, zsep completes but **~174 ZST driver tests fail with ratios ~1e15**.
* With the `sterf` loop fix added, **zsep is identical to reference LAPACK**.
* With the one-token `lqf` fix, **zsvd is identical to reference LAPACK** (42,800 tests, the same single inherent miss).

### 5.4 Examples: numerical validation

* **Across standards:** for every compiler, C++11, C++14 and (GCC) C++17 outputs are bit-identical. In the patched copy, C++11 through C++26 are bit-identical.
* **Across GCC versions:** GCC 14, 15 and 16 are bit-identical.
* **GCC vs. Clang:**
  * 83 programs are identical.
  * 12 complex-arithmetic LAPACK examples differ by at most 1.5e-13 absolute (rounding: libstdc++ vs. libc++ complex arithmetic, FMA contraction).
  * **1 real discrepancy: `tiny_gemv`** (§6.2.2).
* **`-O0` vs. `-O2`:**
  * Identical except `lu_blk-test`, whose residual ‖A−LU‖ is 4.7609e-6 vs. 4.7659e-6 (rounding).
  * `tiny_gemv` under Clang is **correct at `-O0` and NaN at `-O2`**, which is the signature of undefined behaviour.
* **Generic FLENS kernels vs. OpenBLAS CBLAS backend** (`-DWITH_OPENBLAS`; 103 examples + 7 tests):
  * 81 identical, 10 at rounding level.
  * Residuals in `lu_blk`/`lu_unblk` differ at rounding level.
  * `lapack-gels-case2` and `lapack-complex-gels-case2/3` differ substantively. See below.
* **FLENS-LAPACK vs. reference LAPACK** (`-DALWAYS_USE_CXXLAPACK`, OpenBLAS, 56 LAPACK examples):
  * 40 identical, 11 at rounding level (≤ 2e-13).
  * 3 differ: `lapack-gels-case3` and `lapack-complex-gels-case2/3`.
  * **Cause:** these examples call `gels` on the **rank-2** matrix [1 2 3 4; 5 6 7 8; 9 10 11 12] (or its transpose). `lapack-gels.doc` itself states *"It is assumed that A has full rank."* Every backend returns *a* solution of Ax=b, but which one is decided by rounding. FLENS's answer (1,2,3,4) happens to be the true minimum-norm solution (verified: orthogonal to the null space). This is an example-data issue, not a FLENS defect.
* **Apple Accelerate (the documented vecLib path) vs. OpenBLAS:** all 56 LAPACK/tut05 examples identical.
* **Manual spot checks:**
  * `lapack-gesv` gives x = (1,9,9,9); A·x reproduces b exactly.
  * `addingvectors` gives y = 19 everywhere, which is correct.
  * The playground solvers cg, bicgstab and tfqmr agree on (1, −1+i, −2−i).

---

## 6. Detailed failure analysis

Notation: **[C]** compile error, **[L]** link error, **[R]** runtime error/crash, **[N]** numerically wrong result, **[B]** build system, **[D]** external dependency/environment, **[UB]** undefined behaviour.

### 6.1 Blockers for modern C++ standards (FLENS; affect every program)

**6.1.1 [C] `std::allocator<T>::construct`/`destroy` removed in C++20.**
* Locations, 8 sites in 4 files:
  * `flens/storage/fullstorage/fullstorage.tcc:773,784`
  * `flens/storage/packedstorage/packedstorage.tcc:268,281`
  * `flens/storage/array/array.tcc:269,280`
  * `flens/storage/bandstorage/bandstorage.tcc:902,917`
* Configs: all GCC at C++20/23/26; all Clang at C++20 and later.
* Diagnostic: `'std::allocator<double>' has no member named 'construct'` (GCC); `no member named 'destroy' in 'std::allocator<bool>'` (Clang).
* Fix: `std::allocator_traits<decltype(allocator_)>::construct(allocator_, p, value)`, and likewise for `destroy`. This is valid C++11 and verified on all standards.

**6.1.2 [C] `std::ptr_fun` and `std::bind2nd` removed in C++17.**
* Locations, 8 sites:
  * `flens/io/fullstorage/load.tcc:252,320,394,466,540`
  * `flens/io/packedstorage/load.tcc:213,288,361`
* Configs:
  * libc++ (all Clang) at C++17 and later: hard error, `no member named 'bind2nd' in namespace 'std'`.
  * libstdc++ still ships them but warns `-Wdeprecated-declarations` in every mode, including C++11 (16 warnings per translation unit). That is what breaks the `-Werror` interface Makefiles (§6.3.1).
* Fix: `[](char c) { return std::tolower(c, std::locale("")); }`, verified.
* **[Correction]** `std::locale("")` itself throws under most locales with libstdc++ on macOS. The committed fix uses `std::locale::classic()` (A.4).
* Side note: this lower-casing sits inside `#ifndef NDEBUG` and only feeds `ASSERT`s on the MatrixMarket header, so it does not change behaviour.

Applying 6.1.1 and 6.1.2 together gives the "minimal patch". It is 16 replacements across 6 files and produces no output change in C++11.

### 6.2 Genuine library bugs (all compilers, all standards)

**6.2.1 [R][UB] `cxxblas/tinylevel1/copy.tcc:47`: out-of-bounds write.**
* Code: `y[iX] = x[iX];` should be `y[iY] = x[iX];`.
* With `incX=2, incY=1` and n=3, it writes `y[0], y[2], y[4]`, so a 3-element destination is overrun.
* Reproduced: output `1 0 2` instead of `1 2 3`.
* GCC flags `iY` as *set but not used* (`-Wunused-but-set-variable`).

**6.2.2 [N][UB] `cxxblas/tinylevel2/gemv.tcc:50,67,82,101`: reads the result vector when beta = 0.**
* `y[iY] *= beta` is evaluated even for `y = A*x` (beta = 0). `TinyArray`'s default constructor leaves the data uninitialized, so this reads indeterminate values.
* Observed in `examples/tiny_gemv` (A·(1,1,1,1) should be (10,26,30,32)):
  * GCC 14–16 at -O2/-O0: correct.
  * Clang 17/20/23 and Apple Clang at -O2: **`nan nan nan nan`** for both products.
  * Clang at -O0: correct.
* Reference BLAS semantics require that y is not read when beta = 0.
* Also in the same file: `gemv_t` and `gemv_ct` advance `iX += incY` (lines 86 and 105) instead of `incX`. That is latent while strides are 1.
* Fix: if `beta==0`, assign `y = 0` instead of scaling; and use `incX` in those loops.
* Note: the dense generic `cxxblas::gemv` uses `scal_generic(beta, y)`, which propagates NaN/Inf already present in y when beta = 0. That deviates from BLAS, but it is not observable in the examples because FLENS zero-initializes dense arrays.

**6.2.3 [N][R] `flens/lapack/impl/sterf.tcc:172–178`: QL search loop is broken. This is the root cause of wrong Hermitian eigenvalues.**
```cpp
if (l!=lEnd) {
    for (IndexType m=l; l<=lEnd-1; ++l) {   // declares a NEW m, tests/increments l
```
* The reference `DSTERF` has `DO 40 M = L, LEND-1 … M = LEND`. The inner `m` shadows the outer one, and `l` is incremented instead of `m`.
* When `l == lEnd`, the Fortran sets `M = LEND`; FLENS leaves `m` stale.
* **Effect:** `lapack::sterf` returns the input diagonal (sorted) with `info = 0`. Hence `lapack::ev(false, HeMatrix, …)` and Fortran `zheev('N', …)` return **wrong eigenvalues**:
  * d = (2,2,2), e = (1,1) gives 2,2,2 instead of 2∓√2, 2.
  * A random 5×5 Hermitian matrix gives max error 9.09 vs. the `ev(true,…)` path and vs. reference LAPACK.
  * zdrvst test ratios are ~1e15.
* No shipped example calls `ev(false, …)`, which is why the examples did not show it.
* `-Wshadow` reports this line.
* **Fix (verified, zsep becomes identical to reference):** `for (m=l; m<=lEnd-1; ++m) { … }` with the `if (l!=lEnd)` guard removed.

**6.2.4 [N] `flens/lapack/impl/sterf.tcc:88`: `One / SafeMin` uses the enumerator instead of the variable.**
* `SafeMin` is `MachineParameter` = `'S'` = 83, so `safeMax` is 1/83 instead of about 4.5e307.
* As a result `sSaveMax` is about 0.037, so nearly every block gets rescaled. The scaling is undone afterwards, but overflow protection and rounding behaviour change.
* **[Correction]** For blocks with a large norm the downscaling makes small entries subnormal or zero, which gives wrong eigenvalues (A.3).
* Flagged by GCC `-Wdeprecated-enum-float-conversion` (in `-Wall -Wextra` with `-Wdeprecated`). It is a hard error in C++26 with Clang.
* Fix: `safeMin`.

**6.2.5 [N] `flens/lapack/ge/jsv.tcc:907`: `sqrt(small) / Eps` uses the enumerator `Eps`** (`'E'` = 69).
* Reference `DGEJSV` divides by `EPSLN`.
* Wrong perturbation threshold when `l2pert` is set.
* Fix: `eps` (the local declared at line 142).
* Found by inspection plus a grep scan for unshadowed enumerators in arithmetic. Not exercised by the tests run here.
* **[Correction]** The line is inside commented-out code and is never compiled. `jsv`/`dgejsv_` do not work at all (A.3).

**6.2.6 [R] `flens/lapack/interface/src/heev.cc:58`: `zheev_` aborts for N = 0.**
* It builds `DArrayView(3*N-2, RWORK, 1)`, which has length −2 when N=0. The reference LAPACK dimension is `max(1, 3N-2)`.
* Fix (verified): `std::max(INTEGER(1), 3*(*N)-2)`.

**6.2.7 [R] `flens/lapack/ge/lqf.tcc:116`: null pointer passed as the `opts` string.**
* Code: `ilaenv<T>(2, "GELQF", 0, m, n)`. All siblings pass `""`.
* In the interface-library configuration (`USE_ILAENV_WITH_UNDERSCORE`), this reaches `strlen(nullptr)` and crashes `zgelqf_`/`dgelqf_` whenever the workspace is less than optimal. That crashed the zsvd suite.
* Fix (verified, zsvd becomes identical to reference): `""`.

**6.2.8 [R] `flens/lapack/la/lan.tcc:365`: `for (IndexType i=1; i<=min(m,j); ++j)`.**
* The loop increments the outer `j`.
* `lan(InfinityNorm, …)` of an *upper, non-unit* `TrMatrix` therefore sums only row 1 and walks past the last column.
* Reproduced: assertion `col<firstCol_+numCols_` in `fullstorageview.tcc:125`. With `-DNDEBUG` that is an out-of-bounds read.
* Fix: `++i`.

**6.2.9 [R] `flens/lapack/interface/src/orgqr.cc:50` (also `ungqr.cc:49`): workspace query not implemented.**
* `if (lQuery) { // TODO: implement lqf_wsq; ASSERT(0); }`
* Any caller doing the standard `LWORK=-1` query aborts. That is why dcsd/zcsd fail.
* This is a known limitation, not a regression.

**6.2.10 [C] `flens/blas/closures/level1/copy.tcc:160`: `FLENS_BLASLOG_BEGIN_COPY(x1x2, y)`. `x1x2` is undeclared; it should be `x`.**
* Effect: **`-DFLENS_DEBUG_CLOSURES` is unusable with every compiler**, Clang 17 included. That breaks tutorials `tut02-page03-example` and `tut02-page04-example1/2` (and `addingvectors*` in that mode).
* GCC ≥ 14 reports it via `-Wtemplate-body`; Clang via ordinary lookup.
* Fix: `x`.

**6.2.11 [C] `std::complex<int>` produced by index-variable closures (libstdc++ only).**
* `examples/complex-gematrix-views.cc:25` and `trmatrix-indexvariable.cc` write `Z(i,j) = Complex(i+j, j-i)`.
* `flens/scalaroperations/complex.h:50` creates `std::complex<Promotion<int,int>>` = `std::complex<int>`, and assigning that to `complex<double>` fails in libstdc++ (GCC 14, 15 and 16): `cannot convert 'std::complex<int>' to '__complex__ double'`.
* Instantiating `std::complex` with `int` is unspecified behaviour per [complex.numbers]/2.
* libc++ accepts it.
* Fix: promote to the destination's element type, or at least to `double`, in `ScalarOpComplex`.

**6.2.12 [C] `tut08-page01-example2.cc:12`, `tut08-page01-example5.cc:12`: stale API.**
* `TpMatrix<PackedStorage<T, Upper>>` passes a `StorageUpLo` where `PackedStorage` expects a `StorageOrder`. This fails with every compiler.
* There is no `.doc` page for tut08, so this looks like unfinished tutorial code.
* Fix: `TpMatrix<PackedStorage<T> >` and `Matrix Ap(n, Upper);`.

**6.2.13 [C] `flens/hacks/mpfr_real.h:89,101` and `flens/hacks/qd.h:56,68`: user specializations of `std::is_floating_point` / `std::is_arithmetic`.**
* This is UB ([meta.rqmts]).
* Clang ≥ 20 makes it an error (`-Winvalid-specialization`).
* Clang 17 rejects the QD one (`explicit specialization … after instantiation`) and has an ambiguity for MPFR geev.
* GCC: `pow(const mpfr::real&, mpfr::real)` is ambiguous at `cxxblas/auxiliary/pow.tcc:99`.
* Net effect:
  * `mpfr-real-lapack-getrf` works only with GCC (C++11/14/17) and Clang 17.
  * `mpfr-real-lapack-geev` works nowhere.
  * `qd-*` never builds; under GCC that is a link failure, see §6.4.
* Fix: replace the trait hacks with FLENS-owned traits such as `flens::IsReal`, and provide `pow` overloads for these types.

**6.2.14 [C] `playground/flens/solver/cgs.h:33–34`: include-guard typo** (`#ifndef …_CGS_H` / `#define …_CGs_H`).
* `cgs.h` gets included twice, which gives "redefinition of default argument".
* This breaks **all 18 buildable playground programs** with every compiler.
* With this single character fixed (scratch copy), **11/18 build and run** on GCC 16 and Clang 23:
  * bicgstab, cg, cgs, determinant, the fft examples (export/import wisdom, float, matrix, vector), pcg, tfqmr, trace.
* The remaining playground failures are in §6.4 and 6.2.15.

**6.2.15 [C] Other playground and example defects.**
* `example-ptfqmr.cc:56` calls the nonexistent `solver::pbicgstab`.
* `example-suitesparse.cc:72` uses the namespace `suitesparse`, which is never declared.
* `playground/flens/sparse/superlu/superlu.h:38–60` `#include`s SuperLU's C headers *inside* `namespace flens::superlu`. With current libc headers and SuperLU 7 this produces `'memchr' has not been declared in '::'` and similar errors.

**6.2.16 Latent header-guard defects** (no current failure observed):
* Guard and `#define` don't match:
  * `cxxblas/level1extensions/ccopy.h` (`CCOCPY` vs `CCOPY`)
  * `cxxblas/sparselevel3/syccsmm.tcc` (`SYCCSMM` vs `HECCSMM`)
  * `cxxlapack/interface/lae2.tcc`
* `flens/blas/blas.h`, `blas.tcc`, `level1/level1.h` and `level1.tcc` have `#ifndef` but no `#define`.
* **Duplicate guards**, which would silently drop the second header if both were included first-in-TU:
  * `syev.h`/`sygv.h`, `hbgv.h`/`hegv.h`, `syevd.h`/`sygvd.h`, `hbgvd.h`/`hegvd.h`, `lag2d.h`/`lag2s.h`, `la_hercond_c.h`/`la_hercond_x.h`, `her.h`/`her2.h`
  * `heswapr.h`, `hprfs.h`, `la_herpvgrw.h`, `pftrs.h` (all `CXXLAPACK_INTERFACE_DUMMY_H`)
  * `ormtr.h`/`unmtr.h` (and `la_rpvgrw.h`) reuse their `.tcc` guard.

**6.2.17 The 3-argument `lapack::ev(bool, HeMatrix&&, DenseVector&&)` overload does not compile** (`gematrix.h:422: 'std::remove_reference<bool>::type' is not a class`), while the 5-argument form works. Found while writing a reproduction; the root cause was not analysed.

### 6.3 Build-system problems

1. **[B] Hardcoded compiler and `-Werror`.**
   * `flens/blas/interface/blas/Makefile` and `flens/lapack/interface/Makefile` have `CXX = g++-4.8` (exit 127 on any modern system) and `-Wall -Wextra -Werror`.
   * With `CXX=g++-16` or `clang++`, the build fails: 240 / 155 errors (BLAS) and 271 (LAPACK), **all** `-Wdeprecated-declarations` from §6.1.2.
   * With `-Werror` removed it builds cleanly, with only deprecation and `unused-local-typedef` warnings.
2. **[B] Missing targets.**
   * README `make check`: no such target.
   * Top-level `check-lapack` calls `make -C flens/lapack check`: no rule.
   * Top-level `dev-check` calls `make -C flens/blas/interface dev-check`: no rule.
   * The working entry points are `make -C flens/blas/interface check` and `make -C flens/lapack dev-check`.
3. **[B] `flens/lapack/interface/Makefile` `dev` target.** It compiles `*.dev_o` but then runs `ar crs $(FLENSLAPACK) $(OBJECT_FILES)`. So it archives the *non-dev* `.o` files into `lapack_FLENS.a`, and `lapack_FLENS_dev.a` is never produced.
4. **[B] No header dependency tracking.** Verified: after `touch cxxblas/level1/axpy.tcc`, `make -n` rebuilds **0** objects. Incremental builds silently keep stale objects for a header-only library. Touching a `.cc` rebuilds correctly.
5. **[B] Parallel builds.**
   * `make -j4` / `-j8` / `-j10` worked for the reference BLAS, CBLAS and LAPACK builds and for both interface libraries.
   * No races were observed. The `test/` tarball extraction is an order-only directory target.
6. **[B] CMake.**
   * `cmake_minimum_required(VERSION 2.4.8)` is rejected by CMake ≥ 4.0: *Compatibility with CMake < 3.5 has been removed*.
   * With `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` it configures and `cmake --install` installs 2,190 headers. Compiling a tutorial against the installed tree works.
   * Fix: `cmake_minimum_required(VERSION 3.5)`, and move it *before* `project()`.
7. **[B][D] Shipped Netlib Fortran sources vs. modern gfortran.**
   * `dblat1.f` (*Rank mismatch*) and reference LAPACK (*Type mismatch*) fail with gfortran ≥ 10.
   * Building requires `-fallow-argument-mismatch` (passed via `FC=`, `FORTRAN=`, `NOOPT=` on the command line).
8. **[B] Tarball hygiene.**
   * `cxxlapack/netlib/lapack.tgz` contains a prebuilt **x86_64** `xerbla_nu.o`. It is excluded by the Makefile's `filter-out`, so it is harmless here.
   * The top-level `config` file is a stray bare-repository git config.

### 6.4 External-dependency and environment issues (not FLENS defects per se)

* **SuperLU ≥ 5 API:** `examples/gecrs-superlu.cc:89` calls `dgssvx` without the `GlobalLU_t*` parameter added in SuperLU 5.0. The example targets SuperLU 4.3, while Homebrew ships 7.0.1. `geccs-superlu` (which uses `dgssv`) builds and runs fine.
* **QD + GCC:** link errors such as `operator<<(std::ostream&, qd_real const&)` undefined. Homebrew's `libqd` is built against libc++, so it cannot link with libstdc++. That is a stdlib mismatch in the environment. Under Clang the build stops earlier, on §6.2.13.
* **MPI C++ bindings:** `playground/flens/mpi/rank.h:41` uses `MPI::Comm`. The C++ bindings were deprecated in MPI-2.2, removed in MPI-3.0, and are absent from Open MPI 5. Affects all 3 MPI examples.
* **Intel MKL:** not available for macOS/arm64, so `example-mkldss` was skipped.

---

## 7. Warnings, sanitizer findings, numerical discrepancies

### 7.1 Strong warnings

Flags: `-Wall -Wextra -Wpedantic -Wshadow -Wnon-virtual-dtor -Woverloaded-virtual -Wnull-dereference -Wdeprecated`, C++11, all 130 programs. The table shows unique source locations, with the diagnostic count in parentheses.

| Warning | GCC 16 | Clang 23 | Assessment |
|---|---|---|---|
| `-Wshadow` | 208 locs (8,414) | 81 locs (8,024) | mostly benign style; **but it pinpoints the real `sterf.tcc:174` bug**, so worth reviewing in `flens/lapack` |
| `-Wdeprecated-declarations` | 8 (2,048) | 8 (2,000) | `ptr_fun`/`bind2nd` (§6.1.2) |
| `-Wunused-local-typedef(s)` | 9 | 9 | cosmetic |
| `-Wdeprecated-copy-with-user-provided-dtor/-copy` | – | 13 locs | `Matrix`/`Vector`/`Scalar` CRTP bases, tiny views and `ElementClosure` declare destructors or copy-assignment but rely on implicit copy constructors; deprecated since C++11 |
| `-Wdeprecated-enum-float-conversion` | **1** | – | **`sterf.tcc:88` real bug** (§6.2.4) |
| `-Wunused-but-set-variable` | 5 | – | **`tinylevel1/copy.tcc:46` real bug** (§6.2.1); `qp3.tcc:69 iws`, `lsy.tcc:339 lWorkMin`, `examples/permutation.h` benign |
| `-Wdeprecated-variadic-comma-omission` | 1 | – | `flens/auxiliary/isconvertible.h:54`, `(...)` without comma is deprecated in C++26 |
| `-Wimplicit-fallthrough` | 2 | – | `flens/lapack/ge/bal.tcc:127`, intentional (translates Fortran `GOTO`) |
| `-Wdeprecated-dynamic-exception-spec` | – | 3 | `external/real.hpp:61` (`throw(...)` specs, ill-formed in C++17 and later for non-empty specs) |
| `-Wunused-parameter` / `-variable` | 25 / 2 | 30 / 1 | tests and playground, benign; `laqp2.tcc:190 CZero` unused |
| macro redefined | 4 | 4 | `tut05-page02-example2..5` `#define USE_CXXLAPACK` when it is also given with `-D` (harmless) |

### 7.2 Sanitizers

Builds: ASan + UBSan, `-O1 -g`, C++11. GCC 16: 89 clean runs; Clang 23: 69.

* **AddressSanitizer: no reports in any example or test.**
  * The tiny-copy OOB (§6.2.1) is not reached by any shipped program; it was reproduced separately.
  * The uninitialized read in tiny `gemv` (§6.2.2) needs MemorySanitizer, which is not available on macOS.
* **UBSan reports** in 5 runs (GCC) and 26 runs (Clang). All are one pattern, a **biased base pointer that is null for empty views**:
  * `flens/storage/fullstorage/fullstorageview.tcc:134,157`: "applying non-zero offset to null pointer" and "reference binding to address with insufficient space". These come from `operator()` on 0×0 views (21 LAPACK examples, `lu_*`).
  * `cxxblas/sparselevel2/{ge,sy,he}{crs,ccs}mv.tcc:62–69`: `--ia; ja -= ia[1];` on empty matrices. That is the `nnz==0` run of every sparse test.
  * `flens/vectortypes/impl/densevector.tcc:535`: `imag()` of an empty vector.
* FLENS stores pointers pre-offset so it can index from 1. Strictly, any pointer formed before the start of an array is UB, even when non-null. Only the null case is diagnosed.
* No incorrect results were linked to it; all affected runs exit 0 with correct output.
* Risk is **low** but non-zero with aggressive optimizers.

### 7.3 Numerical discrepancies (summary)

| Program / routine | Discrepancy | Cause |
|---|---|---|
| `lapack::ev(false, HeMatrix)` / `zheev('N')` | eigenvalues wrong (error 9.1 on 5×5; test ratios ~1e15) | **bug** §6.2.3 |
| `examples/tiny_gemv` | `nan` vs. correct values (Clang -O2) | **UB** §6.2.2 |
| `lapack-gels-case2/3`, `lapack-complex-gels-case2/3` | backend-dependent solutions | rank-deficient input to `gels` (example data) |
| 12 complex LAPACK examples, GCC vs. Clang | ≤ 1.5e-13 abs | rounding (libstdc++/libc++ complex, FMA) |
| `lu_blk-test`, `lu_unblk-test` residual | ~1e-3 relative on a 4.7e-6 residual | rounding |
| Netlib dsep/dgd/zgd/zsep/zsvd threshold misses | same with pure reference LAPACK | test suite / toolchain, not FLENS |

---

## 8. Recommended fixes, by severity

As originally recommended for `84043d90`. The status column shows what has been done since (Part A).

### Critical: wrong results or memory corruption

| # | Fix | Status |
|---|---|---|
| 1 | `flens/lapack/impl/sterf.tcc:172–178`: QL search loop (§6.2.3) | **fixed** `f90daf7b` |
| 2 | `flens/lapack/impl/sterf.tcc:88`: `SafeMin` → `safeMin` (§6.2.4) | **fixed** `2131b158` |
| 3 | `cxxblas/tinylevel1/copy.tcc:47`: `y[iX]` → `y[iY]` (§6.2.1) | **fixed** `3b9d5be7` (also `ccopy.tcc`) |
| 4 | `cxxblas/tinylevel2/gemv.tcc`: do not read y for beta = 0; `incY` → `incX` (§6.2.2) | **fixed** `5dd06ec5` |
| 5 | `flens/lapack/la/lan.tcc:365`: `++j` → `++i` (§6.2.8) | **fixed** `2da9dc2a` |
| 6 | `flens/lapack/ge/jsv.tcc:907`: `Eps` → `eps` (§6.2.5) | **not changed**: dead code; `jsv`/`dgejsv_` not functional (A.3) |

### High: crashes in the Fortran-callable libraries; modern-standard blockers

| # | Fix | Status |
|---|---|---|
| 7 | `flens/lapack/interface/src/heev.cc:58`: RWORK length `max(1,3N-2)` (§6.2.6) | **fixed** `6d224bd8` |
| 8 | `flens/lapack/ge/lqf.tcc:116`: `0` → `""` (§6.2.7) | **fixed** `98a72f89` |
| 9 | C++20: `std::allocator_traits` (§6.1.1) | **fixed** `9fffd7f6` |
| 10 | C++17/libc++: lambda instead of `bind2nd(ptr_fun(tolower))` (§6.1.2) | **fixed** `9fffd7f6` (with `std::locale::classic()`, `909570db`) |

### Medium: build tooling and broken features

| # | Fix | Status |
|---|---|---|
| 11 | Interface Makefiles: configurable `CXX`, optional `-Werror`, `-MMD -MP`, `dev` archive, `check`/`dev-check` targets | **fixed** `17daf8e5`, `5a710c11` |
| 12 | CMake ≥ 4 | **fixed** `9529905e` |
| 13 | `flens/blas/closures/level1/copy.tcc:160`: `x1x2` → `x` | open |
| 14 | `playground/flens/solver/cgs.h:34`: include guard | open |
| 15 | Netlib Fortran sources: `-fallow-argument-mismatch` | **fixed** `7914b0a1` (probed in `toolchain.mk`, also used for the test suites) |
| 16 | `std::complex<int>` in `ScalarOpComplex` (§6.2.11) | open |

### Low: examples, latent issues, external APIs

All open (17–28 as listed for the baseline):

17. Duplicate or mismatched include guards (§6.2.16).
18. `tut08-page01-example2/5` stale `PackedStorage` usage (§6.2.12).
19. MPFR/QD trait hacks: use FLENS-owned traits (§6.2.13).
20. `orgqr`/`ungqr` workspace query (`LWORK=-1`) (§6.2.9).
21. 3-argument `he::ev` overload (§6.2.17).
22. `gecrs-superlu`: SuperLU ≥ 5 `dgssvx` signature. Playground SuperLU/SuiteSparse/ptfqmr/MPI-C++-bindings updates.
23. Biased-pointer UB on empty views (§7.2).
24. Generic `gemv`/`gemm` with beta = 0 should not propagate NaN from y (BLAS semantics).
25. `lapack-gels-case2/3`: use full-rank example matrices, or point to `gelsy` for the rank-deficient case.
26. `test/geccsmv`: 3000 runs (≈ 20–40 min) where the siblings use 30. Probably a leftover from debugging.
27. Document that `addingvectors*` need `-DUSE_SIMPLE_ET`.
28. Consider enabling `-Wshadow` and `-Wdeprecated-enum-float-conversion` in FLENS-LAPACK development builds; both pointed directly at real bugs.

New items from A.4: MatrixMarket `SyMatrix`/`HeMatrix` loading, packed save/load, `gecopy` ConjTrans, `dgejsv_` LWORK check, the remaining `-Werror` warnings.

---

## 9. Reproducible commands

### 9.1 In the repository (current state)

```sh
make check                                  # regression, Netlib BLAS, Netlib LAPACK (check set)
make check CXX=g++-16 CXXSTD=c++20 FC=gfortran-16
make check-regression REFLAPACK="-L/opt/homebrew/opt/openblas/lib -lopenblas"
make check-lapack-extended                  # all d/z LAPACK tests
make dev-check                              # bitwise comparison with reference LAPACK
make distclean
```

The variables are documented in `toolchain.mk`. Expected deviations of the LAPACK suite are listed in `flens/lapack/check-expected.txt`.

### 9.2 Archive

Everything used for this report is kept in **`/Users/lehn/git/FLENS-compat/`** (`A` below):

| Path | Content |
|---|---|
| `harness/manifest.py` | program inventory and per-program build profiles (§3.2) |
| `harness/run_matrix.py` | parallel compile+run driver → `results.jsonl` |
| `harness/summarize.py`, `final_table.py` | matrix tables |
| `harness/compare_outputs.py` | numeric-tolerant stdout comparison |
| `harness/warnings.py` | warning aggregation |
| `harness/patch_minimal.py` | applies the minimal C++17/20 patch to a copy |
| `harness/run_netlib.sh` | Netlib BLAS + LAPACK (extended set) on an exported tree, writes `summary.txt` |
| `deps/mtl4` | MTL4 checkout for `design-flaws-mtl4` |
| `patches/` | `minimal-cxx17-cxx20.patch`, `lapack-bugfixes.patch` (scratch-copy patches of the investigation) |
| `repro/` | `sterf2.cc`, `heevN.cc`, `tinycopy.cc`, `lantr.cc` |
| `results/baseline-2026-10-26/` | `matrix_O2`, `matrix_patched_O2`, `matrix_O0`, `matrix_san`, `matrix_warn`, `backend`, `matrix_pg`, summaries |
| `results/after-numerical-fixes`, `results/compat-14878d34` | intermediate matrices |
| `results/final-56731834{,-O0,-san,-geccsmv,-regression}` | final validation (A.2) |
| `netlib/baseline-2026-10-26/`, `netlib/rerun-baseline-84043d90/` | Netlib trees and logs of the baseline |
| `netlib/pure-reference-lapack/` | the same suites against pure reference LAPACK |
| `netlib/after-numerical-fixes-98a72f89/`, `netlib/compat-14878d34-{gcc16,clang23}-cxx20/` | after the fixes |
| `netlib/final-56731834-make-targets/` | logs of `make check`, `dev-check`, `check-lapack-extended` and the distclean checks |

### 9.3 Harness runs

The harness does not modify the checkout. Source trees are exported with `git archive`, and `S` is any scratch directory.

```sh
A=/Users/lehn/git/FLENS-compat
git -C /Users/lehn/git/FLENS archive HEAD | tar -x -C $S/src
cd $A/harness
# 7 compilers x C++11..26, -O2 (geccsmv separately, it takes ~20-40 min)
FLENS_ROOT=$S/src python3 -I run_matrix.py --out $S/matrix --opt=-O2 --filter '^(?!test/geccsmv$)'
python3 -I final_table.py $S/matrix/results.jsonl
# -O0, sanitizers, warnings
FLENS_ROOT=$S/src python3 -I run_matrix.py --out $S/O0 --opt=-O0 --compilers gcc16,clang23 --stds 20 --timeout-scale 3
FLENS_ROOT=$S/src python3 -I run_matrix.py --out $S/san --opt="-O1 -g -fno-omit-frame-pointer" \
    --extra="-fsanitize=address,undefined" --compilers gcc16,clang23 --stds 20 --timeout-scale 6
python3 -I run_matrix.py --out $S/warn --extra="-Wall -Wextra -Wpedantic -Wshadow -Wdeprecated" \
    --compilers gcc16,clang23 --stds 11 --no-run
# compare stdout with the baseline
python3 -I compare_outputs.py --only-diff $A/results/baseline-2026-10-26/matrix_O2/gcc16-cxx11 $S/matrix/gcc16-cxx11 $S/matrix/gcc16-cxx20
# backend comparisons (baseline)
OB="-I/opt/homebrew/opt/openblas/include -L/opt/homebrew/opt/openblas/lib -lopenblas"
python3 -I run_matrix.py --out $S/backend --compilers gcc16 --stds 11 --tag=-openblas --extra="-DWITH_OPENBLAS $OB" --filter '^(examples|test/(?!geccsmv$))'
python3 -I run_matrix.py --out $S/backend --compilers gcc16 --stds 11 --tag=-refLAPACK --extra=-DALWAYS_USE_CXXLAPACK --filter '^examples/(lapack-|tut05)'
```

**Baseline by hand** (`84043d90`):
```sh
cd flens/examples
g++-16  -std=c++11 -Wall -I../.. tut01-page01-example.cc -o t && ./t      # works
g++-16  -std=c++20 -I../.. tut01-page01-example.cc -o t                   # baseline: allocator::construct error
clang++ -std=c++17 -I../.. tut01-page01-example.cc -o t                   # baseline: bind2nd error
```

**Bug reproductions** (`A/repro/`; all correct after the fixes):

| File | Bug | Wrong output (baseline) | Expected |
|---|---|---|---|
| `sterf2.cc` | `sterf` (§6.2.3) | `d = 2 2 2` | 0.5858 2 3.4142 |
| `heevN.cc` | `ev(false,…)` (§6.2.3) | max diff 9.09 | ~1e-14 |
| `tinycopy.cc` | tiny `copy` (§6.2.1) | `1 0 2` | `1 2 3` |
| `lantr.cc` | `lan` InfinityNorm (§6.2.8) | assertion | 54 |
| `examples/tiny_gemv` | tiny `gemv` (§6.2.2) | `nan` with `clang++ -O2` | (10, 26, 30, 32) |

---

## 10. What was not tested

* **Other platforms:** only macOS 26 / arm64. No Linux or x86_64 run, so GNU ld group ordering (`LINKER_GROUP_BEGIN`), the `-lstdc++` default in `toolchain.mk` on Linux, and other Linux-specific paths were not exercised.
* **Other compilers:** Intel `icc`/`icx` (listed in the README) were not available. MSVC is not applicable. No GCC older than 14 or Clang older than 17 was used. gfortran 16 was the only Fortran compiler.
* **MemorySanitizer, Valgrind:** not available on macOS/arm64. Uninitialized reads were found only through the GCC/Clang output difference.
* **Netlib LAPACK single/complex-single suites:** not applicable, since `lapack_FLENS.a` exports only d/z routines. The LIN `proto`/RFP/mixed-precision inputs (`dstest`, `zctest`, `*_rfp`) were not run.
* **Netlib suites and compilers:**
  * The extended LAPACK set was run with GCC 16 (C++11, C++20) and Clang 23 (C++20).
  * `check-lapack` was run with Apple Clang 21 (C++11), GCC 16 and Clang 23.
  * `dev-check` was run with GCC 16 (C++11) and Clang 23 (C++20).
  * Other compiler/standard combinations ran the regression tests and the example matrix, but not the Netlib suites.
* **Sanitizers and `-O0`** in the final validation: GCC 16 and Clang 23 at C++20 only. The baseline ran them at C++11.
* **`examples/design-flaws-mtl4`** builds and runs against MTL4 from Git, but it demonstrates an MTL4 behaviour, not FLENS.
* **Playground programs requiring MKL (1) or MPI C++ bindings (3):** could not be built in this environment, so MPI was never run under `mpirun`.
* **`doctool/`** (documentation generator) and the `.doc` → HTML pipeline: not tested.
* **Incremental builds:** checked for the interface libraries and the regression tests; the examples have no build system.
* **CMake:** only configure/build/install with CMake 4.4.4. The install was verified by compiling a tutorial against the installed headers.
