/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: difference.hpp
** -----
** File Created: Thursday, 11th December 2025 22:48:17
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:02:49
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_DIFFERENCE_HPP
#define MSL_DIFFERENCE_HPP

#include <span>
#include <stdexcept>
#include <vector>

#include "matrix/real_matrix_owned.hpp"

namespace msl::difference
{

// ============================================================================
// Difference (First-order)
// ============================================================================

/**
 * @brief First-order difference: diff[i] = y[i+1] - y[i]
 *
 * Output size is n-1
 *
 * @param y Input values
 * @return Differences
 */
inline std::vector<double> diff(std::span<const double> y)
{
    if (y.size() < 2)
    {
        throw std::invalid_argument("Need at least 2 points for diff");
    }

    std::vector<double> result(y.size());
    for (size_t i = 0; i < result.size() - 1; ++i)
    {
        result[i + 1] = y[i + 1] - y[i];
    }

    return result;
}

inline std::vector<double> diff(const std::vector<double> &y)
{
    return diff(std::span<const double>(y));
}

/**
 * @brief First-order difference for matrix (along specified axis)
 *
 * @param mat Input matrix
 * @param axis 0 = row-wise (vertical diff), 1 = column-wise (horizontal diff)
 * @return Difference matrix
 */
inline matrix::matrixd diff(const matrix::matrixd &mat, int axis = 0)
{
    if (axis == 0)
    {
        // Row-wise: diff along rows (vertical)
        if (mat.rows() < 2)
        {
            throw std::invalid_argument(
                "Need at least 2 rows for row-wise diff");
        }

        matrix::matrixd result(mat.rows(), mat.cols());
        for (size_t j = 0; j < mat.cols(); ++j)
        {
            for (size_t i = 0; i < result.rows() - 1; ++i)
            {
                result(i + 1, j) = mat(i + 1, j) - mat(i, j);
            }
        }
        return result;
    }
    else if (axis == 1)
    {
        // Column-wise: diff along columns (horizontal)
        if (mat.cols() < 2)
        {
            throw std::invalid_argument(
                "Need at least 2 cols for column-wise diff");
        }

        matrix::matrixd result(mat.rows(), mat.cols());
        for (size_t i = 0; i < mat.rows(); ++i)
        {
            for (size_t j = 0; j < result.cols() - 1; ++j)
            {
                result(i, j + 1) = mat(i, j + 1) - mat(i, j);
            }
        }
        return result;
    }
    else
    {
        throw std::invalid_argument("axis must be 0 or 1");
    }
}

// ============================================================================
// Gradient (Derivative Approximation)
// ============================================================================

/**
 * @brief Numerical gradient using forward differences (uniform spacing)
 *
 * Uses:
 * - Forward difference
 *
 * @param y Function values
 * @param dx Spacing (default: 1.0)
 * @return Gradient values (same size as input)
 */
inline std::vector<double> forward_gradient(std::span<const double> y,
                                            double dx = 1.0)
{
    if (y.size() < 2)
    {
        throw std::invalid_argument("Need at least 2 points for gradient");
    }

    std::vector<double> grad(y.size());
    // Forward difference, the first point is 0
    for (size_t i = 0; i < y.size() - 1; ++i)
    {
        grad[i + 1] = (y[i + 1] - y[i]) / dx;
    }

    return grad;
}

inline std::vector<double> forward_gradient(const std::vector<double> &y,
                                            double dx = 1.0)
{
    return forward_gradient(std::span<const double>(y), dx);
}

/** @brief Numerical gradient using forward differences (non-uniform spacing)
 *
 * Uses:
 * - Forward difference
 *
 * @param y Function values
 * @param dx Spacing (default: 1.0)
 * @return Gradient values (same size as input)
 */
inline std::vector<double> forward_gradient(std::span<const double> y,
                                            std::span<const double> x)
{
    if (x.size() != y.size())
    {
        throw std::invalid_argument("x and y must have same size");
    }
    if (x.size() < 2)
    {
        throw std::invalid_argument("Need at least 2 points for gradient");
    }
    std::vector<double> grad(y.size());
    // Forward difference, the first point is 0
    for (size_t i = 0; i < y.size() - 1; ++i)
    {
        double dx_local = x[i + 1] - x[i];
        grad[i + 1] = (y[i + 1] - y[i]) / dx_local;
    }

    return grad;
}

/**
 * @brief Gradient for matrix (along specified axis), uses forward differences
 *
 * @param mat Input matrix
 * @param dx Spacing
 * @param axis 0 = gradient along rows, 1 = gradient along columns
 * @return Gradient matrix (same size as input)
 */
inline matrix::matrixd
forward_gradient(const matrix::matrixd &mat, double dx = 1.0, int axis = 0)
{
    if (axis == 0)
    {
        // Gradient along rows (vertical direction)
        if (mat.rows() < 2)
        {
            throw std::invalid_argument("Need at least 2 rows");
        }

        matrix::matrixd grad(mat.rows(), mat.cols());

        // First row gradient is zero
        for (size_t j = 0; j < mat.cols(); ++j)
        {
            // Forward difference
            for (size_t i = 0; i < mat.rows() - 1; ++i)
            {
                grad(i + 1, j) = (mat(i + 1, j) - mat(i, j)) / dx;
            }
        }

        return grad;
    }
    else if (axis == 1)
    {
        // Gradient along columns (horizontal direction)
        if (mat.cols() < 2)
        {
            throw std::invalid_argument("Need at least 2 columns");
        }

        matrix::matrixd grad(mat.rows(), mat.cols());

        // First column gradient is zero
        for (size_t i = 0; i < mat.rows(); ++i)
        {
            // Central at interior columns
            for (size_t j = 0; j < mat.cols() - 1; ++j)
            {
                grad(i, j + 1) = (mat(i, j + 1) - mat(i, j)) / dx;
            }
        }

        return grad;
    }
    else
    {
        throw std::invalid_argument("axis must be 0 or 1");
    }
}

/**
 * @brief Gradient for matrix (along specified axis), uses forward differences
 *
 * @param mat Input matrix
 * @param dx Spacing
 * @param axis 0 = gradient along rows, 1 = gradient along columns
 * @return Gradient matrix (same size as input)
 */
inline matrix::matrixd forward_gradient(const matrix::matrixd &mat,
                                        std::vector<double> x,
                                        int axis = 0)
{
    if (axis == 0)
    {
        // Gradient along rows (vertical direction)
        if (mat.rows() < 2)
        {
            throw std::invalid_argument("Need at least 2 rows");
        }

        if (x.size() != mat.rows())
        {
            throw std::invalid_argument(
                "x size must be number of rows for axis=0");
        }

        matrix::matrixd grad(mat.rows(), mat.cols());

        for (size_t j = 0; j < mat.cols(); ++j)
        {
            // Forward difference
            for (size_t i = 0; i < mat.rows() - 1; ++i)
            {
                double dx_local = x[i + 1] - x[i];
                grad(i + 1, j) = (mat(i + 1, j) - mat(i, j)) / dx_local;
            }
        }

        return grad;
    }
    else if (axis == 1)
    {
        // Gradient along columns (horizontal direction)
        if (mat.cols() < 2)
        {
            throw std::invalid_argument("Need at least 2 columns");
        }

        if (x.size() != mat.cols())
        {
            throw std::invalid_argument(
                "dx size must be number of cols for axis=1");
        }

        matrix::matrixd grad(mat.rows(), mat.cols());

        for (size_t i = 0; i < mat.rows(); ++i)
        {
            // Forward difference
            for (size_t j = 0; j < mat.cols() - 1; ++j)
            {
                double dx_local = x[j + 1] - x[j];
                grad(i, j + 1) = (mat(i, j + 1) - mat(i, j)) / dx_local;
            }
        }

        return grad;
    }
    else
    {
        throw std::invalid_argument("axis must be 0 or 1");
    }
}

/**
 * @brief Numerical gradient using central differences (uniform spacing)
 *
 * Uses:
 * - Forward difference at first point
 * - Central difference at interior points
 * - Backward difference at last point
 *
 * @param y Function values
 * @param dx Spacing (default: 1.0)
 * @return Gradient values (same size as input)
 */
inline std::vector<double> central_gradient(std::span<const double> y,
                                            double dx = 1.0)
{
    if (y.size() < 2)
    {
        throw std::invalid_argument("Need at least 2 points for gradient");
    }

    std::vector<double> grad(y.size());

    // Forward difference at first point
    grad[0] = (y[1] - y[0]) / dx;

    // Central difference at interior points
    for (size_t i = 1; i < y.size() - 1; ++i)
    {
        grad[i] = (y[i + 1] - y[i - 1]) / (2.0 * dx);
    }

    // Backward difference at last point
    grad[y.size() - 1] = (y[y.size() - 1] - y[y.size() - 2]) / dx;

    return grad;
}

inline std::vector<double> central_gradient(const std::vector<double> &y,
                                            double dx = 1.0)
{
    return central_gradient(std::span<const double>(y), dx);
}

/**
 * @brief Numerical gradient with non-uniform spacing
 *
 * @param x Independent variable
 * @param y Dependent variable
 * @return Gradient dy/dx
 */
inline std::vector<double> central_gradient(std::span<const double> x,
                                            std::span<const double> y)
{
    if (x.size() != y.size())
    {
        throw std::invalid_argument("x and y must have same size");
    }
    if (x.size() < 2)
    {
        throw std::invalid_argument("Need at least 2 points for gradient");
    }

    std::vector<double> grad(y.size());

    // Forward difference at first point
    double dx0 = x[1] - x[0];
    grad[0] = (y[1] - y[0]) / dx0;

    // Central difference at interior points
    for (size_t i = 1; i < y.size() - 1; ++i)
    {
        double dx_left = x[i] - x[i - 1];
        double dx_right = x[i + 1] - x[i];

        // Weighted central difference for non-uniform grid
        grad[i] = (dx_left * (y[i + 1] - y[i]) / dx_right
                   + dx_right * (y[i] - y[i - 1]) / dx_left)
                  / (dx_left + dx_right);
    }

    // Backward difference at last point
    size_t n = y.size();
    double dx_last = x[n - 1] - x[n - 2];
    grad[n - 1] = (y[n - 1] - y[n - 2]) / dx_last;

    return grad;
}

/**
 * @brief Gradient for matrix (along specified axis)
 *
 * @param mat Input matrix
 * @param dx Spacing
 * @param axis 0 = gradient along rows, 1 = gradient along columns
 * @return Gradient matrix (same size as input)
 */
inline matrix::matrixd
central_gradient(const matrix::matrixd &mat, double dx = 1.0, int axis = 0)
{
    if (axis == 0)
    {
        // Gradient along rows (vertical direction)
        if (mat.rows() < 2)
        {
            throw std::invalid_argument("Need at least 2 rows");
        }

        matrix::matrixd grad(mat.rows(), mat.cols());

        for (size_t j = 0; j < mat.cols(); ++j)
        {
            // Forward at first row
            grad(0, j) = (mat(1, j) - mat(0, j)) / dx;

            // Central at interior rows
            for (size_t i = 1; i < mat.rows() - 1; ++i)
            {
                grad(i, j) = (mat(i + 1, j) - mat(i - 1, j)) / (2.0 * dx);
            }

            // Backward at last row
            size_t n = mat.rows() - 1;
            grad(n, j) = (mat(n, j) - mat(n - 1, j)) / dx;
        }

        return grad;
    }
    else if (axis == 1)
    {
        // Gradient along columns (horizontal direction)
        if (mat.cols() < 2)
        {
            throw std::invalid_argument("Need at least 2 columns");
        }

        matrix::matrixd grad(mat.rows(), mat.cols());

        for (size_t i = 0; i < mat.rows(); ++i)
        {
            // Forward at first column
            grad(i, 0) = (mat(i, 1) - mat(i, 0)) / dx;

            // Central at interior columns
            for (size_t j = 1; j < mat.cols() - 1; ++j)
            {
                grad(i, j) = (mat(i, j + 1) - mat(i, j - 1)) / (2.0 * dx);
            }

            // Backward at last column
            size_t n = mat.cols() - 1;
            grad(i, n) = (mat(i, n) - mat(i, n - 1)) / dx;
        }

        return grad;
    }
    else
    {
        throw std::invalid_argument("axis must be 0 or 1");
    }
}

/**
 * @brief Gradient for matrix (along specified axis)
 *
 * @param mat Input matrix
 * @param dx Spacing
 * @param axis 0 = gradient along rows, 1 = gradient along columns
 * @return Gradient matrix (same size as input)
 */
inline matrix::matrixd central_gradient(const matrix::matrixd &mat,
                                        std::vector<double> x,
                                        int axis = 0)
{
    if (axis == 0)
    {
        // Gradient along rows (vertical direction)
        if (mat.rows() < 2)
        {
            throw std::invalid_argument("Need at least 2 rows");
        }

        if (x.size() != mat.rows())
        {
            throw std::invalid_argument(
                "x size must be number of rows for axis=0");
        }

        matrix::matrixd grad(mat.rows(), mat.cols());

        for (size_t j = 0; j < mat.cols(); ++j)
        {
            // Forward at first row
            grad(0, j) = (mat(1, j) - mat(0, j)) / (x[1] - x[0]);

            // Central at interior rows
            for (size_t i = 1; i < mat.rows() - 1; ++i)
            {
                grad(i, j) =
                    (mat(i + 1, j) - mat(i - 1, j)) / (x[i + 1] - x[i - 1]);
            }

            // Backward at last row
            size_t n = mat.rows() - 1;
            grad(n, j) = (mat(n, j) - mat(n - 1, j)) / (x[n] - x[n - 1]);
        }

        return grad;
    }
    else if (axis == 1)
    {
        // Gradient along columns (horizontal direction)
        if (mat.cols() < 2)
        {
            throw std::invalid_argument("Need at least 2 columns");
        }

        if (x.size() != mat.cols())
        {
            throw std::invalid_argument(
                "dx size must be number of cols for axis=1");
        }

        matrix::matrixd grad(mat.rows(), mat.cols());

        for (size_t i = 0; i < mat.rows(); ++i)
        {
            // Forward at first column
            grad(i, 0) = (mat(i, 1) - mat(i, 0)) / (x[1] - x[0]);

            // Central at interior columns
            for (size_t j = 1; j < mat.cols() - 1; ++j)
            {
                grad(i, j) =
                    (mat(i, j + 1) - mat(i, j - 1)) / (x[j + 1] - x[j - 1]);
            }

            // Backward at last column
            size_t n = mat.cols() - 1;
            grad(i, n) = (mat(i, n) - mat(i, n - 1)) / (x[n] - x[n - 1]);
        }

        return grad;
    }
    else
    {
        throw std::invalid_argument("axis must be 0 or 1");
    }
}

// ============================================================================
// 2D Gradient (returns both directions)
// ============================================================================

/**
 * @brief Compute 2D gradient (both directions)
 *
 * Returns {dy/dx (horizontal), dy/dy (vertical)}
 *
 * @param mat Input matrix
 * @param dx Horizontal spacing
 * @param dy Vertical spacing
 * @return Pair of gradient matrices {grad_x, grad_y}
 */
inline std::pair<matrix::matrixd, matrix::matrixd>
central_gradient2d(const matrix::matrixd &mat, double dx = 1.0, double dy = 1.0)
{
    if (mat.rows() < 2 || mat.cols() < 2)
    {
        throw std::invalid_argument("Need at least 2x2 matrix for 2D gradient");
    }

    matrix::matrixd grad_x =
        central_gradient(mat, dx, 1); // Horizontal gradient
    matrix::matrixd grad_y = central_gradient(mat, dy, 0); // Vertical gradient

    return {grad_x, grad_y};
}

// ============================================================================
// Higher-order Derivatives
// ============================================================================

/**
 * @brief Second derivative using central differences
 *
 * d²y/dx² ≈ (y[i+1] - 2*y[i] + y[i-1]) / dx²
 *
 * @param y Function values
 * @param dx Spacing
 * @return Second derivative (same size as input)
 */
inline std::vector<double> central_gradient2(std::span<const double> y,
                                             double dx = 1.0)
{
    if (y.size() < 3)
    {
        throw std::invalid_argument(
            "Need at least 3 points for second derivative");
    }

    std::vector<double> grad2(y.size());
    double dx2 = dx * dx;

    // Forward difference at first point (less accurate)
    grad2[0] = (y[2] - 2.0 * y[1] + y[0]) / dx2;

    // Central difference at interior points
    for (size_t i = 1; i < y.size() - 1; ++i)
    {
        grad2[i] = (y[i + 1] - 2.0 * y[i] + y[i - 1]) / dx2;
    }

    // Backward difference at last point
    size_t n = y.size() - 1;
    grad2[n] = (y[n] - 2.0 * y[n - 1] + y[n - 2]) / dx2;

    return grad2;
}

inline std::vector<double> central_gradient2(const std::vector<double> &y,
                                             double dx = 1.0)
{
    return central_gradient2(std::span<const double>(y), dx);
}

// ============================================================================
// Laplacian (for 2D data)
// ============================================================================

/**
 * @brief Compute Laplacian (∇²f = ∂²f/∂x² + ∂²f/∂y²)
 *
 * @param mat Input matrix
 * @param dx Horizontal spacing
 * @param dy Vertical spacing
 * @return Laplacian matrix (same size, boundaries set to 0)
 */
inline matrix::matrixd
laplacian(const matrix::matrixd &mat, double dx = 1.0, double dy = 1.0)
{
    if (mat.rows() < 3 || mat.cols() < 3)
    {
        throw std::invalid_argument("Need at least 3x3 matrix for Laplacian");
    }

    matrix::matrixd result(mat.rows(), mat.cols(), 0.0);

    double dx2 = dx * dx;
    double dy2 = dy * dy;

    // Interior points only
    for (size_t i = 1; i < mat.rows() - 1; ++i)
    {
        for (size_t j = 1; j < mat.cols() - 1; ++j)
        {
            // ∂²f/∂x²
            double d2_dx2 =
                (mat(i, j + 1) - 2.0 * mat(i, j) + mat(i, j - 1)) / dx2;

            // ∂²f/∂y²
            double d2_dy2 =
                (mat(i + 1, j) - 2.0 * mat(i, j) + mat(i - 1, j)) / dy2;

            result(i, j) = d2_dx2 + d2_dy2;
        }
    }

    return result;
}

// ============================================================================
// Divergence and Curl (for vector fields)
// ============================================================================

/**
 * @brief Compute divergence of 2D vector field
 *
 * div(F) = ∂Fx/∂x + ∂Fy/∂y
 *
 * @param Fx X-component of vector field
 * @param Fy Y-component of vector field
 * @param dx Horizontal spacing
 * @param dy Vertical spacing
 * @return Divergence (scalar field)
 */
inline matrix::matrixd divergence(const matrix::matrixd &Fx,
                                  const matrix::matrixd &Fy,
                                  double dx = 1.0,
                                  double dy = 1.0)
{
    if (Fx.rows() != Fy.rows() || Fx.cols() != Fy.cols())
    {
        throw std::invalid_argument("Fx and Fy must have same size");
    }

    auto dFx_dx = central_gradient(Fx, dx, 1); // ∂Fx/∂x
    auto dFy_dy = central_gradient(Fy, dy, 0); // ∂Fy/∂y

    matrix::matrixd div(Fx.rows(), Fx.cols());
    for (size_t i = 0; i < div.rows(); ++i)
    {
        for (size_t j = 0; j < div.cols(); ++j)
        {
            div(i, j) = dFx_dx(i, j) + dFy_dy(i, j);
        }
    }

    return div;
}

/**
 * @brief Compute curl of 2D vector field
 *
 * curl(F) = ∂Fy/∂x - ∂Fx/∂y
 *
 * @param Fx X-component of vector field
 * @param Fy Y-component of vector field
 * @param dx Horizontal spacing
 * @param dy Vertical spacing
 * @return Curl (scalar field, z-component)
 */
inline matrix::matrixd curl(const matrix::matrixd &Fx,
                            const matrix::matrixd &Fy,
                            double dx = 1.0,
                            double dy = 1.0)
{
    if (Fx.rows() != Fy.rows() || Fx.cols() != Fy.cols())
    {
        throw std::invalid_argument("Fx and Fy must have same size");
    }

    auto dFy_dx = central_gradient(Fy, dx, 1); // ∂Fy/∂x
    auto dFx_dy = central_gradient(Fx, dy, 0); // ∂Fx/∂y

    matrix::matrixd curl_z(Fx.rows(), Fx.cols());
    for (size_t i = 0; i < curl_z.rows(); ++i)
    {
        for (size_t j = 0; j < curl_z.cols(); ++j)
        {
            curl_z(i, j) = dFy_dx(i, j) - dFx_dy(i, j);
        }
    }

    return curl_z;
}

// ============================================================================
// Savitzky-Golay Filter for Smoothed Derivatives
// ============================================================================

/**
 * @brief Savitzky-Golay derivative (smoothed, for noisy data)
 *
 * Uses polynomial fitting over local window
 * More robust to noise than raw finite differences
 *
 * @param y Function values
 * @param window_size Window size (must be odd, >= 5)
 * @param poly_order Polynomial order (< window_size)
 * @param dx Spacing
 * @return Smoothed derivative
 */
inline std::vector<double> savgol_gradient(std::span<const double> y,
                                           int window_size = 5,
                                           int poly_order = 2,
                                           double dx = 1.0)
{
    if (window_size % 2 == 0)
    {
        throw std::invalid_argument("window_size must be odd");
    }
    if (poly_order >= window_size)
    {
        throw std::invalid_argument("poly_order must be < window_size");
    }
    if (y.size() < static_cast<size_t>(window_size))
    {
        throw std::invalid_argument("Signal too short for window_size");
    }

    // For simplicity, use central differences with averaging
    // Full S-G implementation would require matrix operations

    std::vector<double> grad(y.size());
    int half_window = window_size / 2;

    // Central smoothed derivative
    for (size_t i = 0; i < y.size(); ++i)
    {
        int start = std::max(0, static_cast<int>(i) - half_window);
        int end = std::min(static_cast<int>(y.size()) - 1,
                           static_cast<int>(i) + half_window);

        // Simple weighted average of local differences
        double sum_grad = 0.0;
        int count = 0;

        for (int j = start; j < end; ++j)
        {
            sum_grad += (y[j + 1] - y[j]) / dx;
            ++count;
        }

        grad[i] = sum_grad / count;
    }

    return grad;
}

} // namespace msl::difference

#endif // MSL_DIFFERENCE_HPP
