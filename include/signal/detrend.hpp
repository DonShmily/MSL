/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: detrend.hpp
** -----
** File Created: Sunday, 14th December 2025 17:39:55
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:40:49
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_DETREND_HPP
#define MSL_DETREND_HPP

#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

#include "polynomial/polynomial.hpp"

namespace msl::signal
{
// Detrend data by removing polynomial trend of degree n
inline std::vector<double> detrend(std::span<const double> data,
                                   std::size_t n = 1)
{
    if (data.size() <= n)
    {
        throw std::invalid_argument(
            "Data size must be greater than polynomial degree");
    }

    // Generate x values
    std::vector<double> x(data.size());
    std::iota(x.begin(), x.end(), 0.0);

    // Fit polynomial to data
    msl::polynomial::Polynomial poly(x, data, n);

    // Subtract trend
    std::vector<double> detrended(data.size());
    for (std::size_t i = 0; i < data.size(); ++i)
    {
        detrended[i] = data[i] - poly(x[i]);
    }

    return detrended;
}

// Detrend data by removing polynomial trend of degree n
inline std::vector<double> detrend(std::span<const double> x,
                                   std::span<const double> data,
                                   std::size_t n = 1)
{
    if (x.size() != data.size())
    {
        throw std::invalid_argument("x and data must have same size");
    }

    if (data.size() <= n)
    {
        throw std::invalid_argument(
            "Data size must be greater than polynomial degree for detrending");
    }

    // Fit polynomial to data
    msl::polynomial::Polynomial poly(x, data, n);

    // Subtract trend
    std::vector<double> detrended(data.size());
    for (std::size_t i = 0; i < data.size(); ++i)
    {
        detrended[i] = data[i] - poly(x[i]);
    }

    return detrended;
}

// Detrend data by removing polynomial trend of degree n
inline matrix::matrixd detrend(const matrix::matrixd &data, std::size_t n = 1)
{
    matrix::matrixd detrended(data.rows(), data.cols());
    for (std::size_t col = 0; col < data.cols(); ++col)
    {
        // Extract column data
        auto col_data = data.column(col);

        // Detrend column
        auto detrended_col = detrend(col_data, n);

        // Store detrended column
        std::copy(detrended_col.begin(),
                  detrended_col.end(),
                  detrended.data() + col * data.rows());
    }
    return detrended;
}

} // namespace msl::signal

#endif // MSL_DETREND_HPP