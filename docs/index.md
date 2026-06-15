# MSL — Modern Scientific Library

## Overview

MSL is a **C++20 header-only** numerical computing library inspired by MATLAB and NumPy/SciPy. It provides a comprehensive set of tools for matrix operations, signal processing, numerical differentiation, integration, interpolation, polynomial fitting, scalar nonlinear equations, and ODE initial value problems — all without external runtime dependencies beyond Eigen3.

### Design Philosophy

- **Header-only**: Zero compilation overhead for distribution. Include what you need.
- **Column-major storage**: Fortran/MATLAB-compatible memory layout (`data_[j * rows_ + i]`).
- **Zero-copy API**: Every computational function provides two overloads — one that writes to a caller-provided buffer, and one that allocates and returns a new vector.
- **Eigen3 bridge**: Seamless zero-copy interop with Eigen for high-performance linear algebra and FFT.
- **No virtual dispatch overhead**: Matrix hierarchy uses CRTP-like patterns with `std::span` rather than virtual methods.

### Modules

| Module | Namespace | Aggregator | Purpose |
|--------|-----------|------------|---------|
| [matrix](matrix.md) | `msl::matrix` | `msl/matrix.hpp` | Real/complex matrices, views, linear algebra, decompositions |
| [signal](signal.md) | `msl::signal` | `msl/signal.hpp` | FFT, filter design, Butterworth IIR, zero-phase filtering, PSD, window functions |
| [difference](difference.md) | `msl::difference` | `msl/difference.hpp` | Diff, gradient, Laplacian, divergence, curl, Savitzky-Golay |
| [equation](equation.md) | `msl::equation` | `msl/equation.hpp` | One-dimensional nonlinear equation solvers |
| [integral](integral.md) | `msl::integral` | `msl/integral.hpp` | Trapezoidal, Simpson, Romberg, adaptive quadrature |
| [interp](interp.md) | `msl::interp` | `msl/interp.hpp` | 1D interpolation (linear, spline, Akima, PCHIP, nearest, polynomial) |
| [ode](ode.md) | `msl::ode` | `msl/ode.hpp` | Fixed-step and adaptive ODE initial value solvers |
| [polynomial](polynomial.md) | `msl::polynomial` | `msl/polynomial.hpp` | Polynomial class, least-squares fitting, evaluation |
| [utils](utils.md) | `msl::utils` | — | File I/O, assert handler |

## Quick Start

```cpp
#include "msl/matrix.hpp"
#include "msl/difference.hpp"
#include "msl/integral.hpp"
#include "msl/interp.hpp"
#include "msl/equation.hpp"
#include "msl/ode.hpp"
#include "msl/signal/fft.hpp"
#include "msl/signal/filter.hpp"
#include "msl/signal/filtfilt.hpp"

// Create a 3x3 matrix (row-major initializer list → column-major storage)
msl::matrix::matrixd A(3, 3, {1, 2, 3, 4, 5, 6, 7, 8, 9});

// Transpose
auto At = msl::matrix::transpose(A);

// Matrix multiplication
auto C = A * At;

// QR decomposition
auto [Q, R] = msl::matrix::qr(A);

// Difference
std::vector<double> y = {0.0, 1.0, 4.0, 9.0, 16.0};
auto dy = msl::difference::diff(y); // {1.0, 3.0, 5.0, 7.0}

// Integral
auto integral = msl::integral::trapz(y, 1.0); // ≈ 21.5

// Interpolation
msl::interp::Linear interp;
interp.set_data({0.0, 1.0, 2.0}, {0.0, 1.0, 4.0});
double val = interp(1.5); // ≈ 2.5

// Polynomial fitting
auto coeffs = msl::polynomial::polyfit(
    std::vector{0.0, 1.0, 2.0, 3.0},
    std::vector{0.0, 1.0, 4.0, 9.0}, 2);
// coeffs ≈ [0, 0, 1] → f(x) = x^2

// Scalar nonlinear equation
double root = msl::equation::brent_root(
    [](double x) { return x * x - 2.0; }, 0.0, 2.0);

// ODE initial value problem: y' = y, y(0) = 1
auto y1 = msl::ode::rk4_final(
    [](double, const msl::ode::state& y) { return msl::ode::state{y[0]}; },
    std::vector<double>{1.0}, 0.0, 1.0, 0.01);
```

## Dependency

Only **Eigen3** is required. On Linux, install via `libeigen3-dev`. On Windows, use vcpkg manifest mode (`vcpkg.json`).

## License

MSL is distributed under the MIT license. Eigen3, the sole dependency, is licensed under MPL 2.0. See `licenses/LICENSE` and `NOTICE.md`.
