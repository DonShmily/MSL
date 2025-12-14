/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: romberg.hpp
** -----
** File Created: Wednesday, 15th October 2025 22:18:38
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:03:14
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_ROMBERG_HPP
#define MSL_ROMBERG_HPP

#include <cmath>
#include <span>
#include <stdexcept>
#include <vector>

#include "trapz.hpp"

namespace msl::integral
{

// ============================================================================
// Total Integral (Romberg Integration)
// ============================================================================

/**
 * @brief Romberg integration (adaptive, high accuracy)
 *
 * Uses Richardson extrapolation on trapezoidal rule
 *
 * @param y Function values (size must be 2^k + 1)
 * @param dx Uniform spacing
 * @param tol Tolerance for convergence
 * @return Integral value
 */
inline double romberg(std::span<const double> y, double dx, double tol = 1e-10)
{
    size_t n = y.size();

    // Check if n = 2^k + 1
    size_t k = 0;
    size_t check = 1;
    while (check < n)
    {
        check *= 2;
        ++k;
    }
    if (check + 1 != n)
    {
        throw std::invalid_argument("Romberg needs 2^k + 1 points");
    }

    // Romberg table
    std::vector<std::vector<double>> R(k + 1);
    for (size_t i = 0; i <= k; ++i)
    {
        R[i].resize(i + 1);
    }

    // R[0,0] = trapezoidal with all points
    R[0][0] = trapz(y, dx);

    // Fill Romberg table
    size_t stride = 1;
    for (size_t i = 1; i <= k; ++i)
    {
        stride *= 2;

        // Trapezoidal with stride
        double sum = 0.5 * (y.front() + y.back());
        for (size_t j = stride; j < n; j += stride)
        {
            sum += y[j];
        }
        R[i][0] = sum * dx * stride;

        // Richardson extrapolation
        double power = 4.0;
        for (size_t j = 1; j <= i; ++j)
        {
            R[i][j] = (power * R[i][j - 1] - R[i - 1][j - 1]) / (power - 1.0);
            power *= 4.0;
        }

        // Check convergence
        if (i > 0 && std::abs(R[i][i] - R[i - 1][i - 1]) < tol)
        {
            return R[i][i];
        }
    }

    return R[k][k];
}
} // namespace msl::integral

#endif // MSL_ROMBERG_HPP