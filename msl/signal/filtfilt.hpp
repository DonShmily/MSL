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

#ifndef MSL_FILTFILT_HPP
#define MSL_FILTFILT_HPP

#include <algorithm>
#include <span>
#include <stdexcept>
#include <vector>

#include "filter_design.hpp"

namespace msl::signal
{

// ============================================================================
// Forward-backward filtering (zero-phase)
// ============================================================================

inline std::vector<std::vector<double>>
matrix_inverse(const std::vector<std::vector<double>> &mat)
{
    size_t n = mat.size();
    std::vector<std::vector<double>> A = mat;
    std::vector<std::vector<double>> inv(n, std::vector<double>(n, 0.0));

    // Initialize inv as identity matrix
    for (size_t i = 0; i < n; ++i)
    {
        inv[i][i] = 1.0;
    }

    // Gaussian elimination
    for (size_t i = 0; i < n; ++i)
    {
        // Find pivot
        size_t max_row = i;
        double max_val = std::abs(A[i][i]);
        for (size_t k = i + 1; k < n; ++k)
        {
            if (std::abs(A[k][i]) > max_val)
            {
                max_val = std::abs(A[k][i]);
                max_row = k;
            }
        }

        if (max_row != i)
        {
            std::swap(A[i], A[max_row]);
            std::swap(inv[i], inv[max_row]);
        }

        // Normalize current row
        double pivot = A[i][i];
        for (size_t j = 0; j < n; ++j)
        {
            A[i][j] /= pivot;
            inv[i][j] /= pivot;
        }

        // Eliminate
        for (size_t k = 0; k < n; ++k)
        {
            if (k != i)
            {
                double factor = A[k][i];
                for (size_t j = 0; j < n; ++j)
                {
                    A[k][j] -= factor * A[i][j];
                    inv[k][j] -= factor * inv[i][j];
                }
            }
        }
    }

    return inv;
}

inline std::vector<double>
matrix_vector_multiply(const std::vector<std::vector<double>> &mat,
                       const std::vector<double> &vec)
{
    size_t rows = mat.size();
    size_t cols = vec.size();
    std::vector<double> result(rows, 0.0);

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            result[i] += mat[i][j] * vec[j];
        }
    }

    return result;
}

/**
 * @brief Compute initial conditions for filtfilt
 *
 * Computes initial filter state to minimize startup transients
 * Based on MATLAB's filtfilt implementation
 */
inline std::vector<double> compute_filtfilt_zi(const FilterCoefficients &coeffs)
{
    size_t nfilt = std::max(coeffs.a.size(), coeffs.b.size());
    size_t n = nfilt - 1;

    if (n == 0)
    {
        return {};
    }

    //  Extend coefficient vectors
    std::vector<double> a = coeffs.a;
    std::vector<double> b = coeffs.b;
    a.resize(nfilt, 0.0);
    b.resize(nfilt, 0.0);

    // Build sparse matrix sp
    std::vector<int> rows, cols;
    std::vector<double> data;

    // rows = [0:n-1, 1:n-1, 0:n-2]
    for (size_t i = 0; i < n; ++i)
        rows.push_back(i);
    if (n > 1)
    {
        for (size_t i = 1; i < n; ++i)
            rows.push_back(i);
        for (size_t i = 0; i < n - 1; ++i)
            rows.push_back(i);
    }

    // cols = [0, 0, ..., 0, 1:n-1, 1:n-1]
    for (size_t i = 0; i < n; ++i)
        cols.push_back(0);
    if (n > 1)
    {
        for (size_t i = 1; i < n; ++i)
            cols.push_back(i);
        for (size_t i = 1; i < n; ++i)
            cols.push_back(i);
    }

    // data = [1+a[1], a[2:n], ones(1,n-1), -ones(1,n-1)]
    data.push_back(1.0 + a[1]);
    for (size_t i = 2; i < nfilt; ++i)
    {
        data.push_back(a[i]);
    }
    if (n > 1)
    {
        for (size_t i = 0; i < n - 1; ++i)
            data.push_back(1.0);
        for (size_t i = 0; i < n - 1; ++i)
            data.push_back(-1.0);
    }

    //  Build full matrix
    std::vector<std::vector<double>> sp(n, std::vector<double>(n, 0.0));
    for (size_t k = 0; k < rows.size(); ++k)
    {
        sp[rows[k]][cols[k]] += data[k];
    }

    // Compute right-hand side vector: b[1:] - b[0] * a[1:]
    std::vector<double> rhs(n);
    for (size_t i = 0; i < n; ++i)
    {
        rhs[i] = b[i + 1] - b[0] * a[i + 1];
    }

    // Solve: zi = inv(sp) * rhs
    auto sp_inv = matrix_inverse(sp);
    return matrix_vector_multiply(sp_inv, rhs);
}

/**
 * @brief Apply filter with initial conditions
 */
inline std::vector<double> filter_with_zi(std::span<const double> signal,
                                          const FilterCoefficients &coeffs,
                                          std::vector<double> zi)
{
    if (coeffs.a.empty())
    {
        throw std::invalid_argument("Feedback filter coefficients are empty");
    }

    // Normalize coefficients (ensure a[0] = 1.0)
    std::vector<double> a = coeffs.a;
    std::vector<double> b = coeffs.b;

    double a0 = a[0];
    if (a0 == 0.0)
    {
        throw std::invalid_argument(
            "First feedback coefficient must be non-zero");
    }

    if (a0 != 1.0)
    {
        for (auto &val : a)
            val /= a0;
        for (auto &val : b)
            val /= a0;
    }

    const size_t input_size = signal.size();
    size_t filter_order = std::max(a.size(), b.size());

    a.resize(filter_order, 0.0);
    b.resize(filter_order, 0.0);
    zi.resize(filter_order, 0.0);

    std::vector<double> output(input_size);

    // Filtering implementation
    for (size_t i = 0; i < input_size; ++i)
    {
        size_t order = filter_order - 1;
        while (order > 0)
        {
            if (i >= order)
            {
                zi[order - 1] = b[order] * signal[i - order]
                                - a[order] * output[i - order] + zi[order];
            }
            --order;
        }
        output[i] = b[0] * signal[i] + zi[0];
    }

    return output;
}

/**
 * @brief Apply filter forward and backward (zero-phase filtering)
 *
 * This doubles the filter order and eliminates phase distortion.
 * Equivalent to MATLAB's filtfilt() function.
 *
 * Implementation based on MATLAB R2023a filtfilt.m
 * Uses proper initial conditions and edge handling
 *
 * @param signal Input signal
 * @param coeffs Filter coefficients
 * @return Zero-phase filtered signal
 *
 * @note Effective filter order is doubled
 * @note Uses reflection padding with initial conditions
 */
inline std::vector<double> filtfilt(std::span<const double> signal,
                                    const FilterCoefficients &coeffs)
{
    const int len = static_cast<int>(signal.size());
    const int nfilt =
        static_cast<int>(std::max(coeffs.b.size(), coeffs.a.size()));
    const int nfact = 3 * (nfilt - 1);

    if (len <= nfact)
    {
        throw std::invalid_argument(
            "Input data too short! Must have length > 3 * filter_order");
    }

    // Compute initial conditions
    auto zi_base = compute_filtfilt_zi(coeffs);

    //  Left padding: 2*signal[0] - signal[nfact:1:-1]
    std::vector<double> leftpad;
    for (int i = nfact; i >= 1; --i)
    {
        leftpad.push_back(2.0 * signal[0] - signal[i]);
    }

    //  Right padding: 2*signal[end] - signal[end-2:end-nfact-1:-1]
    std::vector<double> rightpad;
    for (int i = len - 2; i >= len - nfact - 1; --i)
    {
        rightpad.push_back(2.0 * signal[len - 1] - signal[i]);
    }

    // Construct padded signal
    std::vector<double> signal1;
    signal1.reserve(leftpad.size() + len + rightpad.size());
    signal1.insert(signal1.end(), leftpad.begin(), leftpad.end());
    signal1.insert(signal1.end(), signal.begin(), signal.end());
    signal1.insert(signal1.end(), rightpad.begin(), rightpad.end());

    // Forward filtering
    std::vector<double> zi = zi_base;
    double y0 = signal1[0];
    for (auto &z : zi)
    {
        z *= y0;
    }
    auto signal2 = filter_with_zi(signal1, coeffs, zi);

    // Reverse
    std::reverse(signal2.begin(), signal2.end());

    // Backward filtering
    zi = zi_base;
    y0 = signal2[0];
    for (auto &z : zi)
    {
        z *= y0;
    }
    signal1 = filter_with_zi(signal2, coeffs, zi);

    // Reverse back
    std::vector<double> result;
    result.reserve(len);
    for (int i = signal1.size() - nfact - 1; i >= nfact; --i)
    {
        result.push_back(signal1[i]);
    }

    return result;
}

} // namespace msl::signal

#endif // MSL_FILTFILT_HPP