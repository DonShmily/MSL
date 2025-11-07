/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025
**  File: window.hpp
**  Description: Apply window functions to signals
*/

#ifndef MSL_WINDOW_HPP
#define MSL_WINDOW_HPP

#include <cmath>
#include <span>
#include <stdexcept>
#include <vector>

#include "matrix/real_matrix_owned.hpp"

namespace msl::signal
{
// ========================================================================
// 1. Generate window functions
// ========================================================================
// Hamming window
inline std::vector<double> hamming_window(size_t n)
{
    std::vector<double> window(n);

    for (size_t i = 0; i < n; ++i)
    {
        window[i] =
            0.54 - 0.46 * std::cos(2.0 * std::numbers::pi * i / (n - 1));
    }

    return window;
}
// Hann window
inline std::vector<double> hann_window(size_t n)
{
    std::vector<double> window(n);

    for (size_t i = 0; i < n; ++i)
    {
        window[i] =
            0.5 * (1.0 - std::cos(2.0 * std::numbers::pi * i / (n - 1)));
    }

    return window;
}
// Blackman window
inline std::vector<double> blackman_window(size_t n)
{
    std::vector<double> window(n);

    for (size_t i = 0; i < n; ++i)
    {
        window[i] = 0.42 - 0.5 * std::cos(2.0 * std::numbers::pi * i / (n - 1))
                    + 0.08 * std::cos(4.0 * std::numbers::pi * i / (n - 1));
    }

    return window;
}
// Rectangular window
inline std::vector<double> rectangular_window(size_t n)
{
    return std::vector<double>(n, 1.0);
}

// ========================================================================
// 2. Apply window functions
// ========================================================================
// Apply window to 1-D signal
inline std::vector<double> apply_window(std::span<const double> signal,
                                        const std::vector<double> &window)
{
    if (signal.size() != window.size())
    {
        throw std::invalid_argument(
            "Signal and window must have the same length");
    }

    std::vector<double> windowed_signal(signal.size());

    for (size_t i = 0; i < signal.size(); ++i)
    {
        windowed_signal[i] = signal[i] * window[i];
    }

    return windowed_signal;
}
// Apply window to 2-D matrix (column-wise)
inline matrixd apply_window_columns(const matrixd &input,
                                    const std::vector<double> &window)
{
    size_t n_rows = input.rows();
    size_t n_cols = input.cols();

    if (n_rows != window.size())
    {
        throw std::invalid_argument(
            "Number of rows in input must match window length");
    }

    matrixd output(n_rows, n_cols);

    for (size_t j = 0; j < n_cols; ++j)
    {
        for (size_t i = 0; i < n_rows; ++i)
        {
            output(i, j) = input(i, j) * window[i];
        }
    }

    return output;
}

} // namespace msl::signal

#endif // MSL_WINDOW_HPP