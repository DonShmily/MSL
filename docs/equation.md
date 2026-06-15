# Equation Module

**Namespace**: `msl::equation`  
**Aggregator**: `msl/equation.hpp`

## Overview

The equation module solves one-dimensional nonlinear equations of the form
`f(x) = 0`. Algorithms return a diagnostic `root_result` by default and also
provide `*_root` convenience wrappers for code that only needs the root value.

## Result Type

```cpp
enum class root_status {
    converged,
    max_iterations,
    invalid_interval,
    zero_derivative,
    non_finite_value
};

struct root_result {
    double root;
    double residual;
    size_t iterations;
    root_status status;

    bool converged() const;
    explicit operator bool() const;
    double value() const; // throws if not converged
};
```

Use the result-returning API when diagnostics matter:

```cpp
auto r = msl::equation::brent(
    [](double x) { return x * x - 2.0; }, 0.0, 2.0);

if (r) {
    double x = r.root;
}
```

Use the convenience API for direct mathematical-style calls:

```cpp
double x = msl::equation::brent_root(
    [](double x) { return x * x - 2.0; }, 0.0, 2.0);
```

## Options

```cpp
struct root_options {
    double tol = 1e-12;
    size_t max_iter = 100;
    double derivative_tol = 1e-14;
};
```

Every algorithm supports either explicit tolerance arguments or `root_options`.

## Algorithms

### Bisection

Robust bracketing method. Requires `f(a)` and `f(b)` to have opposite signs.

```cpp
auto r = msl::equation::bisection(f, a, b);
auto r = msl::equation::bisection(f, a, b, options);
double x = msl::equation::bisection_root(f, a, b);
```

### Brent

Hybrid bracketing method combining bisection, secant, and inverse quadratic
interpolation. This is the recommended default scalar root finder when a
bracket is available.

```cpp
auto r = msl::equation::brent(f, a, b);
auto r = msl::equation::brent(f, a, b, options);
double x = msl::equation::brent_root(f, a, b);
```

### Newton

Fast local method requiring an explicit derivative.

```cpp
auto r = msl::equation::newton(f, df, x0);
auto r = msl::equation::newton(f, df, x0, options);
double x = msl::equation::newton_root(f, df, x0);
```

### Secant

Derivative-free local method using two initial guesses.

```cpp
auto r = msl::equation::secant(f, x0, x1);
auto r = msl::equation::secant(f, x0, x1, options);
double x = msl::equation::secant_root(f, x0, x1);
```

### Regula Falsi

False-position bracketing method. Requires a valid sign-changing interval.

```cpp
auto r = msl::equation::regula_falsi(f, a, b);
auto r = msl::equation::regula_falsi(f, a, b, options);
double x = msl::equation::regula_falsi_root(f, a, b);
```

## MATLAB Validation

Run the C++ test first to generate comparison data:

```sh
xmake run test_equation
```

Then run:

```sh
matlab -batch "run('matlab/validation/equation/validate_equation.m')"
```

The script compares MSL roots against MATLAB `fzero` and known analytic values.
