/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \include\integral\cumtrapz.hpp
** -----
** File Created: Wednesday, 15th October 2025 22:26:59
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Wednesday, 15th October 2025 22:27:28
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_CUMTRAPZ_HPP
#define MSL_CUMTRAPZ_HPP

#include <span>
#include <stdexcept>
#include <vector>

#include "matrix/real_matrix_owned.hpp"

namespace msl::integral
{

// ============================================================================
// Cumulative Trapezoidal Integration
// ============================================================================

/**
 * @brief Cumulative trapezoidal integration for vector (uniform spacing)
 *
 * Computes cumulative integral: output[i] = ∫[0 to i] y dx
 *
 * @param y Function values at equally spaced points
 * @param dx Spacing between points
 * @return Cumulative integral values
 */
inline std::vector<double> cumtrapz(std::span<const double> y, double dx)
{
    if (y.size() < 2)
    {
        throw std::invalid_argument("Need at least 2 points for integration");
    }

    std::vector<double> result(y.size(), 0.0);

    for (size_t i = 1; i < y.size(); ++i)
    {
        result[i] = result[i - 1] + 0.5 * (y[i] + y[i - 1]) * dx;
    }

    return result;
}

inline std::vector<double> cumtrapz(const std::vector<double> &y, double dx)
{
    return cumtrapz(std::span<const double>(y), dx);
}

/**
 * @brief Cumulative trapezoidal integration for vector (non-uniform spacing)
 *
 * @param x Independent variable values
 * @param y Function values at x points
 * @return Cumulative integral values
 */
inline std::vector<double> cumtrapz(std::span<const double> x,
                                    std::span<const double> y)
{
    if (x.size() != y.size())
    {
        throw std::invalid_argument("x and y must have same size");
    }
    if (x.size() < 2)
    {
        throw std::invalid_argument("Need at least 2 points for integration");
    }

    std::vector<double> result(y.size(), 0.0);

    for (size_t i = 1; i < y.size(); ++i)
    {
        double dx = x[i] - x[i - 1];
        result[i] = result[i - 1] + 0.5 * (y[i] + y[i - 1]) * dx;
    }

    return result;
}

/**
 * @brief Cumulative trapezoidal integration for matrix (column-wise, uniform
 * spacing)
 *
 * Integrates each column independently
 *
 * @param mat Input matrix (each column is a function)
 * @param dx Spacing between rows
 * @return Matrix of cumulative integrals
 */
inline matrix::matrixd cumtrapz(const matrix::matrixd &mat, double dx)
{
    if (mat.rows() < 2)
    {
        throw std::invalid_argument("Need at least 2 rows for integration");
    }

    matrix::matrixd result(mat.rows(), mat.cols(), 0.0);

    for (size_t j = 0; j < mat.cols(); ++j)
    {
        for (size_t i = 1; i < mat.rows(); ++i)
        {
            result(i, j) =
                result(i - 1, j) + 0.5 * (mat(i, j) + mat(i - 1, j)) * dx;
        }
    }

    return result;
}

/**
 * @brief Cumulative trapezoidal integration for matrix (column-wise,
 * non-uniform)
 *
 * @param x Independent variable values (same for all columns)
 * @param mat Input matrix
 * @return Matrix of cumulative integrals
 */
inline matrix::matrixd cumtrapz(std::span<const double> x,
                                const matrix::matrixd &mat)
{
    if (x.size() != mat.rows())
    {
        throw std::invalid_argument("x size must match number of rows");
    }
    if (mat.rows() < 2)
    {
        throw std::invalid_argument("Need at least 2 rows for integration");
    }

    matrix::matrixd result(mat.rows(), mat.cols(), 0.0);

    for (size_t j = 0; j < mat.cols(); ++j)
    {
        for (size_t i = 1; i < mat.rows(); ++i)
        {
            double dx = x[i] - x[i - 1];
            result(i, j) =
                result(i - 1, j) + 0.5 * (mat(i, j) + mat(i - 1, j)) * dx;
        }
    }

    return result;
}

} // namespace msl::integral

#endif // MSL_CUMTRAPZ_HPP