# MATLAB Validation

This folder contains MATLAB cross-validation scripts for MSL modules.

The workflow is:

1. Run the corresponding C++ test to generate data under `test_result/<module>/matlab_compare/`.
2. Run the MATLAB validation script.

Examples:

```sh
xmake run test_equation
matlab -batch "run('matlab/validation/equation/validate_equation.m')"

xmake run test_ode
matlab -batch "run('matlab/validation/ode/validate_ode.m')"
```

Generated comparison data is intentionally kept under `test_result/`, which is
ignored by git.

Available validation scripts:

- `difference/validate_difference.m`
- `equation/validate_equation.m`
- `integral/validate_integral.m`
- `interp/validate_interp.m`
- `matrix/validate_matrix.m`
- `ode/validate_ode.m`
- `polynomial/validate_polynomial.m`
- `signal/validate_signal.m`
