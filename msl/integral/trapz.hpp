/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: trapz.hpp
** -----
** File Created: Friday, 9th January 2026 14:58:10
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Thursday, 5th March 2026 14:57:47
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_TRAPZ_HPP
#define MSL_TRAPZ_HPP

#include <span>
#include <stdexcept>
#include <vector>

#include "matrix/real_matrix_base.hpp"

namespace msl::integral
{

// ============================================================================
// Total Integral (Trapezoidal Rule)
// ============================================================================

/**
 * @brief Compute total integral using trapezoidal rule (uniform spacing)
 *
 * @param y Function values
 * @param dx Spacing
 * @return Total integral value
 */
inline double trapz(std::span<const double> y, double dx)
{
    if (y.size() < 2)
    {
        throw std::invalid_argument(
            "Trapz: need at least 2 points for integration");
    }

    double sum = 0.5 * (y.front() + y.back());
    for (size_t i = 1; i < y.size() - 1; ++i)
    {
        sum += y[i];
    }

    return sum * dx;
}

/**
 * @brief Compute total integral using trapezoidal rule (non-uniform spacing)
 *
 * @param x Independent variable values
 * @param y Function values at x points
 * @return Total integral value
 */
inline double trapz(std::span<const double> x, std::span<const double> y)
{
    if (x.size() != y.size())
    {
        throw std::invalid_argument("Trapz: x and y must have same size");
    }
    if (x.size() < 2)
    {
        throw std::invalid_argument(
            "Trapz: need at least 2 points for integration");
    }

    double sum = 0.0;
    for (size_t i = 1; i < y.size(); ++i)
    {
        sum += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
    }

    return sum;
}

/**
 * @brief Compute total integral using trapezoidal rule (uniform spacing)
 *
 * Convenience wrapper that allocates and returns a vector.
 * For zero-copy operations, use the void version with span parameter.
 *
 * @param y Function values
 * @param dx Spacing
 * @return Total integral value
 */
inline void
trapz(const matrix::real_matrix_base &mat, std::span<double> result, double dx)
{
    if (mat.rows() < 2)
    {
        throw std::invalid_argument(
            "Trapz: need at least 2 rows for integration");
    }
    if (result.size() != mat.cols())
    {
        throw std::invalid_argument(
            "Trapz: result span must have same size as number of columns");
    }

    for (size_t j = 0; j < mat.cols(); ++j)
    {
        double sum = 0.5 * (mat(0, j) + mat(mat.rows() - 1, j));
        for (size_t i = 1; i < mat.rows() - 1; ++i)
        {
            sum += mat(i, j);
        }
        result[j] = sum * dx;
    }
}

/**
 * @brief Compute total integral using trapezoidal rule (uniform spacing)
 *
 * Convenience wrapper that allocates and returns a vector. For zero-copy
 * operations, use the void version with span parameter.
 *
 * @param mat Input matrix (each column is a function)
 * @param dx Spacing between rows
 * @return Vector of integral values for each column
 */
inline std::vector<double> trapz(const matrix::real_matrix_base &mat, double dx)
{
    std::vector<double> result(mat.cols());
    trapz(mat, result, dx);
    return result;
}

} // namespace msl::integral

#endif // MSL_TRAPZ_HPP
