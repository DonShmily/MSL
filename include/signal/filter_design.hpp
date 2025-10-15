/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025
**  File: filter_design.hpp
**  Description: Digital filter design (IIR filters)
*/

#ifndef MSL_FILTER_DESIGN_HPP
#define MSL_FILTER_DESIGN_HPP

#include <vector>

namespace msl::signal
{

// ============================================================================
// Filter specifications
// ============================================================================

/**
 * @brief Filter type based on frequency response
 */
enum class FilterType
{
    lowpass,  // Low-pass filter
    highpass, // High-pass filter
    bandpass, // Band-pass filter
    bandstop  // Band-stop (notch) filter
};

/**
 * @brief Filter design method
 */
enum class FilterMethod
{
    butterworth, // Maximally flat passband
    chebyshev1,  // Equiripple passband, monotonic stopband (future)
    chebyshev2,  // Monotonic passband, equiripple stopband (future)
    elliptic     // Equiripple in both bands (future)
};

/**
 * @brief Filter coefficients (numerator and denominator)
 *
 * Represents digital filter as transfer function:
 * H(z) = (b[0] + b[1]*z^-1 + ... + b[M]*z^-M) / (a[0] + a[1]*z^-1 + ... +
 * a[N]*z^-N)
 *
 * Where a[0] is typically normalized to 1.0
 */
struct FilterCoefficients
{
    std::vector<double> b; // Numerator coefficients
    std::vector<double> a; // Denominator coefficients (a[0] = 1.0)

    FilterCoefficients() = default;
    FilterCoefficients(std::vector<double> num, std::vector<double> den)
        : b(std::move(num)), a(std::move(den))
    {}

    [[nodiscard]] size_t order() const
    {
        return std::max(a.size(), b.size()) - 1;
    }
};

} // namespace msl::signal

#endif // MSL_FILTER_DESIGN_HPP