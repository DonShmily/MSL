/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: adaptive_simpson_quadrature.hpp
** -----
** File Created: Friday, 9th January 2026 14:58:10
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Thursday, 5th March 2026 14:58:32
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_ADAPTIVE_SIMPSON_QUADRATURE_HPP
#define MSL_ADAPTIVE_SIMPSON_QUADRATURE_HPP

#include <cmath>
#include <functional>
#include <stdexcept>

namespace msl::integral {
// ============================================================================
// Adaptive Integration (Adaptive Simpson's)
// ============================================================================

namespace detail {
inline double adaptive_simpson_recursive(const std::function<double(double)> &f,
                                         double a,
                                         double b,
                                         double fa,
                                         double fb,
                                         double fc,
                                         double S_ab,
                                         double tol,
                                         int max_depth,
                                         int depth = 0) {
    if (depth >= max_depth) {
        return S_ab;
    }

    double c = 0.5 * (a + b);
    double h = 0.25 * (b - a);

    double fd = f(a + h);
    double fe = f(b - h);

    // Simpson's rule on [a,c] and [c,b]
    double S_ac = h / 3.0 * (fa + 4.0 * fd + fc);
    double S_cb = h / 3.0 * (fc + 4.0 * fe + fb);
    double S_acb = S_ac + S_cb;

    // Error estimate
    double error = (S_acb - S_ab) / 15.0;

    if (std::abs(error) < tol) {
        return S_acb + error;
    }

    // Recursive subdivision
    return adaptive_simpson_recursive(
               f, a, c, fa, fd, fc, S_ac, tol / 2, max_depth, depth + 1)
           + adaptive_simpson_recursive(
               f, c, b, fc, fe, fb, S_cb, tol / 2, max_depth, depth + 1);
}
} // namespace detail

/**
 * @brief Adaptive Simpson's quadrature for function
 *
 * @param f Function to integrate
 * @param a Lower limit
 * @param b Upper limit
 * @param tol Tolerance
 * @param max_depth Maximum recursion depth
 * @return Integral value
 */
inline double quad(const std::function<double(double)> &f,
                   double a,
                   double b,
                   double tol = 1e-8,
                   int max_depth = 50) {
    if (a >= b) {
        throw std::invalid_argument(
            "Quad: Lower limit must be less than upper limit.");
    }

    if (tol <= 0) {
        throw std::invalid_argument("Quad: Tolerance must be positive.");
    }

    if (max_depth <= 0) {
        throw std::invalid_argument("Quad: Max depth must be positive.");
    }

    if (f == nullptr) {
        throw std::invalid_argument("Quad: Function cannot be null.");
    }

    double fa = f(a);
    double fb = f(b);
    double fc = f(0.5 * (a + b));

    double h = b - a;
    double S = h / 6.0 * (fa + 4.0 * fc + fb);

    return detail::adaptive_simpson_recursive(
        f, a, b, fa, fb, fc, S, tol, max_depth);
}

} // namespace msl::integral

#endif // MSL_ADAPTIVE_SIMPSON_QUADRATURE_HPP