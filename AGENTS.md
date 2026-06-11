# MSL — Agent Guide

## Build

Two build systems coexist:

- **xmake** (cross-platform, default platform: mingw):
  `xmake f -p mingw && xmake` — builds only test executables, not library targets. SDK path (`C:/Programing/msys64/ucrt64`) is hardcoded in root `xmake.lua`; adjust for your MSYS2 install.
- **Visual Studio 2022** (Windows only): open `MSL.sln`. Uses v143 toolchain, C++20 (`stdcpp20`).

Single external dependency: `eigen3` via vcpkg manifest mode.

## Run Tests

```sh
xmake run test_matrix        # or test_signal, test_difference, test_integral, test_interp, test_polynomial
```

- No test framework — custom macros: `TEST_CASE`, `EXPECT_TRUE`, `EXPECT_EQ`, `EXPECT_NEAR`, `EXPECT_CPLX_NEAR`.
- Tests write output to `test_result/<module>/` (gitignored).
- Prefer deterministic unit assertions over output-file-only checks. Legacy signal/interp/integral/difference tests may read sample files, but new coverage should use small in-memory fixtures where practical.

## Structure

| Path | Purpose |
|---|---|
| `msl/<module>.hpp` | Aggregator headers (public API entrypoints) |
| `msl/<module>/` | Header-only implementations (C++20, `inline`) |
| `project/` | Visual Studio `.vcxproj` static lib targets |
| `tests/` | xmake test executables |
| `matlab/` | MATLAB cross-validation scripts |

Namespaces mirror directories: `msl::matrix`, `msl::signal`, `msl::difference`, `msl::integral`, `msl::interp`, `msl::polynomial`, `msl::utils`.

**Column-major** storage (`data_[j * rows_ + i]`).

## Gotchas

- **Typo filename**: `msl/matrix/martrix_decompose.hpp` (not `matrix`).
- **Signal aggregator incomplete**: `msl/signal.hpp` does NOT include `filtfilt.hpp` or `filter.hpp` — include them explicitly if needed.
- **C++20 for builds, C++23 for indexing**: xmake uses `c++20`; clangd and VS Code IntelliSense use `c++23`/`gnu++23`.
- **Eigen bridge**: `msl/matrix/eigen_interface.hpp` provides zero-copy `Eigen::Map` views and `matmul()`, `matvec()`, `solve()` wrappers. FFT uses `unsupported/Eigen/FFT`.
- **Format**: `.clang-format` (LLVM-based, column 80, indent 4). **Lint**: `.clang-tidy` (bugprone, cppcoreguidelines, modernize, performance, readability).
- **Header template**: psi-header with author "Dong Feiyue" `<FeiyueDong@outlook.com>`.
- **License**: MIT (MSL) + MPL 2.0 (Eigen3 dependency). See `licenses/LICENSE`.
