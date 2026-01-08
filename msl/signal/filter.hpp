/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: filter_apply.hpp
** -----
** File Created: Tuesday, 14th October 2025 11:12:27
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:05:23
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_FILTER_HPP
#define MSL_FILTER_HPP

#include <algorithm>
#include <span>
#include <stdexcept>
#include <vector>

#include "filter_design.hpp"

namespace msl::signal
{

// ============================================================================
// Direct Form II Transposed (most numerically stable)
// ============================================================================

/**
 * @brief Apply IIR filter to signal using Direct Form II Transposed
 *
 * This is the most numerically stable form for IIR filters.
 *
 * @param signal Input signal
 * @param coeffs Filter coefficients (a[0] must be 1.0)
 * @return Filtered signal
 *
 * @note Uses zero initial conditions
 */
inline std::vector<double> filter(std::span<const double> signal,
                                  const FilterCoefficients &coeffs)
{
    if (coeffs.a.empty() || coeffs.b.empty())
    {
        throw std::invalid_argument("Filter coefficients cannot be empty");
    }

    if (std::abs(coeffs.a[0] - 1.0) > 1e-10)
    {
        throw std::invalid_argument(
            "First denominator coefficient must be 1.0");
    }

    const size_t n = signal.size();
    const size_t na = coeffs.a.size();
    const size_t nb = coeffs.b.size();
    const size_t nz = std::max(na, nb) - 1; // Number of delay states

    std::vector<double> output(n);
    std::vector<double> z(nz, 0.0); // State vector (delay line)

    // Direct Form II Transposed implementation
    for (size_t i = 0; i < n; ++i)
    {
        double x = signal[i];

        // Output = b[0]*x + z[0]
        output[i] = (nb > 0 ? coeffs.b[0] : 0.0) * x + (nz > 0 ? z[0] : 0.0);

        // Update state vector
        for (size_t j = 0; j < nz - 1; ++j)
        {
            double bj = (j + 1 < nb) ? coeffs.b[j + 1] : 0.0;
            double aj = (j + 1 < na) ? coeffs.a[j + 1] : 0.0;
            z[j] = bj * x - aj * output[i] + z[j + 1];
        }

        // Last state
        if (nz > 0)
        {
            double bn = (nz < nb) ? coeffs.b[nz] : 0.0;
            double an = (nz < na) ? coeffs.a[nz] : 0.0;
            z[nz - 1] = bn * x - an * output[i];
        }
    }

    return output;
}

} // namespace msl::signal

#endif // MSL_FILTER_APPLY_HPP