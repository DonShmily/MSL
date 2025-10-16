/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \include\integral\trapz.hpp
** -----
** File Created: Wednesday, 15th October 2025 22:32:51
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Wednesday, 15th October 2025 22:33:01
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_TRAPZ_HPP
#define MSL_TRAPZ_HPP

#include <span>
#include <stdexcept>
#include <vector>

#include "matrix/real_matrix_owned.hpp"

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
        throw std::invalid_argument("Need at least 2 points for integration");
    }

    double sum = 0.5 * (y.front() + y.back());
    for (size_t i = 1; i < y.size() - 1; ++i)
    {
        sum += y[i];
    }

    return sum * dx;
}

inline double trapz(const std::vector<double> &y, double dx)
{
    return trapz(std::span<const double>(y), dx);
}

/**
 * @brief Compute total integral using trapezoidal rule (non-uniform spacing)
 */
inline double trapz(std::span<const double> x, std::span<const double> y)
{
    if (x.size() != y.size())
    {
        throw std::invalid_argument("x and y must have same size");
    }
    if (x.size() < 2)
    {
        throw std::invalid_argument("Need at least 2 points for integration");
    }

    double sum = 0.0;
    for (size_t i = 1; i < y.size(); ++i)
    {
        sum += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
    }

    return sum;
}

/**
 * @brief Integrate each column of matrix and return vector of results
 */
inline std::vector<double> trapz(const matrixd &mat, double dx)
{
    if (mat.rows() < 2)
    {
        throw std::invalid_argument("Need at least 2 rows for integration");
    }

    std::vector<double> result(mat.cols());

    for (size_t j = 0; j < mat.cols(); ++j)
    {
        double sum = 0.5 * (mat(0, j) + mat(mat.rows() - 1, j));
        for (size_t i = 1; i < mat.rows() - 1; ++i)
        {
            sum += mat(i, j);
        }
        result[j] = sum * dx;
    }

    return result;
}

} // namespace msl::integral

#endif // MSL_TRAPZ_HPP
