# Developer Guide

## Build Systems

### xmake (Cross-Platform)

```sh
# Configure for your platform
xmake f -p linux       # Linux (GCC)
xmake f -p mingw       # MinGW (default)
xmake f -p windows     # Windows (MSVC)
xmake f -p macosx      # macOS (Clang)

# Build all test executables
xmake

# Run specific tests
xmake run test_matrix
xmake run test_signal
xmake run test_difference
xmake run test_integral
xmake run test_interp
xmake run test_polynomial
```

- Linux requires `libeigen3-dev` system package
- MinGW has a hardcoded SDK path in `xmake.lua` (`C:/Programing/msys64/ucrt64`)

### Visual Studio 2022 (Windows Only)

Open `MSL.sln`. Uses v143 toolchain with C++20 (`stdcpp20`). Project files are in `project/` and `project/test/`. Both static library targets and test projects are configured.

## Testing

### Custom Test Framework

MSL uses a minimal custom test framework with macros defined in each test file:

```cpp
TEST_CASE("Test name") {
    // test body
    EXPECT_TRUE(condition);
    EXPECT_EQ(expected, actual);
    EXPECT_NEAR(expected, actual, tolerance);
    EXPECT_CPLX_NEAR(expected, actual, tolerance);
}
```

### Running Tests

```sh
# Via xmake
xmake run test_matrix

# Direct execution
./out/linux/bin/test_matrix
```

Tests write output to `test_result/<module>/` (gitignored). The test `main()` function returns the count of failed tests (0 = all passed).

### Test Structure

Each test module (`tests/test_<module>/`) contains:
- `test_<module>.cpp` — main test file
- `xmake.lua` — xmake build config for the test target

The main `tests/xmake.lua` includes all subdirectories.

## Code Style

### Formatting

Clang-format is configured in `.clang-format`:
- Based on LLVM style
- Indent width: 4
- Column limit: 80
- Access modifier offset: -4
- Pointer alignment: right

```sh
find msl tests -name '*.hpp' -o -name '*.cpp' | xargs clang-format -i
```

### Linting

Clang-tidy checks are configured in `.clang-tidy`:
- Enabled: `bugprone-*`, `cppcoreguidelines-*`, `modernize-*`, `performance-*`, `readability-*`
- Specific checks are disabled (`avoid-magic-numbers`, `macro-usage`, etc.)

```sh
find msl tests -name '*.hpp' -o -name '*.cpp' | xargs clang-tidy --fix
```

### Conventions

1. **Header-only**: All implementations marked `inline`, no `.cpp` files in `msl/`
2. **Zero-copy API**: Computational functions provide two overloads:
   - `func(input, output_span)` — caller provides output buffer
   - `func(input)` — convenience, returns newly allocated container
3. **Error handling**: `std::invalid_argument` for input validation, `std::out_of_range` for bounds errors, `assert()` for internal invariants
4. **Doxygen comments**: `/** @brief ... @param ... @return ... @throws ... */`
5. **Column-major storage**: `data_[j * rows_ + i]`
6. **psi-header template**: Author "Dong Feiyue" `<FeiyueDong@outlook.com>`

## Directory Structure

```
msl/
├── <module>.hpp           # Aggregator header (public API)
├── <module>/              # Implementation headers
├── .clang-format          # Format rules
├── .clang-tidy            # Lint rules
├── xmake.lua             # Root xmake config
├── vcpkg.json            # Windows dependencies
├── MSL.sln              # Visual Studio solution
├── project/              # VS project files
├── tests/                # xmake test targets
│   ├── xmake.lua        # Tests include
│   ├── test_<module>/   # Per-module tests
│   └── test_<module>.cpp
├── matlab/               # MATLAB validation scripts
├── resource/             # MATLAB coefficient generation
├── docs/                 # Module documentation
└── licenses/             # MPL 2.0 license text
```

## Adding a New Module

1. Create `msl/<module>.hpp` aggregator that includes all implementation headers
2. Create `msl/<module>/` directory with implementation headers
3. Add namespace `msl::<module>`
4. Create test file `tests/test_<module>/test_<module>.cpp`
5. Add test target `tests/test_<module>/xmake.lua`
6. Update `tests/xmake.lua` to include the new test directory
7. For VS builds, create project files in `project/<module>/` and `project/test/test_<module>/`

## Known Issues

1. **Typo filename**: `msl/matrix/martrix_decompose.hpp` should be `matrix_decompose.hpp`
2. **Incomplete signal aggregator**: `msl/signal.hpp` does not include `filtfilt.hpp` or `filter.hpp`
3. **C++ standard mismatch**: xmake builds with `c++20`; clangd and VS Code IntelliSense use `c++23`/`gnu++23`
