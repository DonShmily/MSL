/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025
**  File: filter_apply.hpp
**  Description: Apply digital filters to signals
*/

#ifndef MSL_FILTER_APPLY_HPP
#define MSL_FILTER_APPLY_HPP

#include <algorithm>
#include <span>
#include <stdexcept>
#include <vector>
#include "filter_design.hpp"
#include "matrix/real_matrix_owned.hpp"


namespace msl::filter
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

/**
 * @brief Apply filter to signal (vector overload)
 */
inline std::vector<double> filter(const std::vector<double> &signal,
                                  const FilterCoefficients &coeffs)
{
    return filter(std::span<const double>(signal), coeffs);
}

// ============================================================================
// Forward-backward filtering (zero-phase)
// ============================================================================

/**
 * @brief Apply filter forward and backward (zero-phase filtering)
 *
 * This doubles the filter order and eliminates phase distortion.
 * Equivalent to MATLAB's filtfilt() function.
 *
 * @param signal Input signal
 * @param coeffs Filter coefficients
 * @return Zero-phase filtered signal
 *
 * @note Effective filter order is doubled
 * @note Uses reflection padding at boundaries
 */
inline std::vector<double> filtfilt(std::span<const double> signal,
                                    const FilterCoefficients &coeffs)
{
    if (signal.size() < 3 * (coeffs.a.size() - 1))
    {
        throw std::invalid_argument(
            "Signal too short for forward-backward filtering");
    }

    // Pad signal with reflections
    const size_t pad_len = 3 * (coeffs.a.size() - 1);
    std::vector<double> padded(signal.size() + 2 * pad_len);

    // Reflect at start
    for (size_t i = 0; i < pad_len; ++i)
    {
        padded[i] = 2.0 * signal[0] - signal[pad_len - i];
    }

    // Copy signal
    std::copy(signal.begin(), signal.end(), padded.begin() + pad_len);

    // Reflect at end
    for (size_t i = 0; i < pad_len; ++i)
    {
        size_t idx = signal.size() - 1 - i;
        padded[pad_len + signal.size() + i] = 2.0 * signal.back() - signal[idx];
    }

    // Forward filtering
    auto forward = filter(padded, coeffs);

    // Reverse
    std::reverse(forward.begin(), forward.end());

    // Backward filtering
    auto backward = filter(forward, coeffs);

    // Reverse again
    std::reverse(backward.begin(), backward.end());

    // Remove padding
    std::vector<double> result(signal.size());
    std::copy(backward.begin() + pad_len,
              backward.begin() + pad_len + signal.size(),
              result.begin());

    return result;
}

/**
 * @brief Zero-phase filtering (vector overload)
 */
inline std::vector<double> filtfilt(const std::vector<double> &signal,
                                    const FilterCoefficients &coeffs)
{
    return filtfilt(std::span<const double>(signal), coeffs);
}

// ============================================================================
// Matrix filtering (column-wise)
// ============================================================================

/**
 * @brief Apply filter to each column of a matrix
 *
 * @param signals Matrix where each column is a signal
 * @param coeffs Filter coefficients
 * @return Filtered matrix
 */
inline matrixd filter_columns(const matrixd &signals,
                              const FilterCoefficients &coeffs)
{
    matrixd output(signals.rows(), signals.cols());

    for (size_t j = 0; j < signals.cols(); ++j)
    {
        auto col_span = signals.column(j);
        auto filtered = filter(col_span, coeffs);

        for (size_t i = 0; i < signals.rows(); ++i)
        {
            output(i, j) = filtered[i];
        }
    }

    return output;
}

/**
 * @brief Apply zero-phase filter to each column of a matrix
 */
inline matrixd filtfilt_columns(const matrixd &signals,
                                const FilterCoefficients &coeffs)
{
    matrixd output(signals.rows(), signals.cols());

    for (size_t j = 0; j < signals.cols(); ++j)
    {
        auto col_span = signals.column(j);
        auto filtered = filtfilt(col_span, coeffs);

        for (size_t i = 0; i < signals.rows(); ++i)
        {
            output(i, j) = filtered[i];
        }
    }

    return output;
}

// ============================================================================
// Convenience functions
// ============================================================================

/**
 * @brief Design and apply lowpass filter in one step
 *
 * @param signal Input signal
 * @param order Filter order
 * @param fc Cutoff frequency (normalized)
 * @param zero_phase Use zero-phase filtering (default: true)
 * @return Filtered signal
 */
inline std::vector<double> lowpass(std::span<const double> signal,
                                   int order,
                                   double fc,
                                   bool zero_phase = true)
{
    auto coeffs = butterworth_lowpass(order, fc);
    return zero_phase ? filtfilt(signal, coeffs) : filter(signal, coeffs);
}

/**
 * @brief Design and apply highpass filter
 */
inline std::vector<double> highpass(std::span<const double> signal,
                                    int order,
                                    double fc,
                                    bool zero_phase = true)
{
    auto coeffs = butterworth_highpass(order, fc);
    return zero_phase ? filtfilt(signal, coeffs) : filter(signal, coeffs);
}

/**
 * @brief Design and apply bandpass filter
 */
inline std::vector<double> bandpass(std::span<const double> signal,
                                    int order,
                                    double fc_low,
                                    double fc_high,
                                    bool zero_phase = true)
{
    auto coeffs = butterworth_bandpass(order, fc_low, fc_high);
    return zero_phase ? filtfilt(signal, coeffs) : filter(signal, coeffs);
}

// Vector overloads
inline std::vector<double> lowpass(const std::vector<double> &signal,
                                   int order,
                                   double fc,
                                   bool zero_phase = true)
{
    return lowpass(std::span<const double>(signal), order, fc, zero_phase);
}

inline std::vector<double> highpass(const std::vector<double> &signal,
                                    int order,
                                    double fc,
                                    bool zero_phase = true)
{
    return highpass(std::span<const double>(signal), order, fc, zero_phase);
}

inline std::vector<double> bandpass(const std::vector<double> &signal,
                                    int order,
                                    double fc_low,
                                    double fc_high,
                                    bool zero_phase = true)
{
    return bandpass(
        std::span<const double>(signal), order, fc_low, fc_high, zero_phase);
}

} // namespace msl::filter

#endif // MSL_FILTER_APPLY_HPP