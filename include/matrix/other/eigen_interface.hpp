/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025
**  File: eigen_interface.hpp
**  Description: Interface between MSL matrix types and Eigen
*/

#ifndef MSL_EIGEN_INTERFACE_HPP
#define MSL_EIGEN_INTERFACE_HPP

#include <eigen3/Eigen/Dense>
#include "matrix/matrix_real.hpp"
#include "matrix/matrix_view.hpp"

namespace msl::eigen_interface
{

// ============================================================================
// 1. To Eigen (from MSL)
// ============================================================================

/// @brief Convert matrix_real (owning) to Eigen::Map
/// @note No data is copied, this is a lightweight view
inline Eigen::Map<Eigen::MatrixXd> as_eigen(matrix_real &mat)
{
    return Eigen::Map<Eigen::MatrixXd>(mat.data(), mat.rows(), mat.cols());
}

/// @brief Convert const matrix_real (owning) to Eigen::Map
inline Eigen::Map<const Eigen::MatrixXd> as_eigen(const matrix_real &mat)
{
    return Eigen::Map<const Eigen::MatrixXd>(
        mat.data(), mat.rows(), mat.cols());
}

/// @brief Convert matrix_view (non-owning) to Eigen::Map
inline Eigen::Map<Eigen::MatrixXd> as_eigen(matrix_view<double> view)
{
    return Eigen::Map<Eigen::MatrixXd>(view.data(), view.rows(), view.cols());
}

/// @brief Convert const matrix_view (non-owning) to Eigen::Map
inline Eigen::Map<const Eigen::MatrixXd>
as_eigen(const matrix_view<double> &view)
{
    return Eigen::Map<const Eigen::MatrixXd>(
        view.data(), view.rows(), view.cols());
}

// ============================================================================
// 2. From Eigen (to MSL)
// ============================================================================

/// @brief Create a new matrix_real from Eigen::MatrixXd (deep copy)
inline matrix_real from_eigen(const Eigen::MatrixXd &eig)
{
    matrix_real result(eig.rows(), eig.cols());
    std::copy(eig.data(), eig.data() + eig.size(), result.data());
    return result;
}

/// @brief Create a matrix_view (non-owning) from Eigen::MatrixXd
/// @warning The view is only valid while Eigen matrix still exists!
inline matrix_view<double> view_from_eigen(Eigen::MatrixXd &eig)
{
    return matrix_view<double>(eig.data(), eig.rows(), eig.cols());
}

/// @brief Const version (read-only)
inline matrix_view<const double> view_from_eigen(const Eigen::MatrixXd &eig)
{
    return matrix_view<const double>(eig.data(), eig.rows(), eig.cols());
}

// ============================================================================
// 3. Convenience helpers
// ============================================================================

/// @brief Copy from Eigen to existing MSL matrix (resize if needed)
inline void copy_from_eigen(matrix_real &dest, const Eigen::MatrixXd &eig)
{
    if (dest.rows() != static_cast<size_t>(eig.rows())
        || dest.cols() != static_cast<size_t>(eig.cols()))
    {
        dest.resize(eig.rows(), eig.cols());
    }
    std::copy(eig.data(), eig.data() + eig.size(), dest.data());
}

/// @brief Copy from MSL matrix to Eigen matrix (resize if needed)
inline void copy_to_eigen(const matrix_real &src, Eigen::MatrixXd &eig)
{
    eig.resize(src.rows(), src.cols());
    std::copy(src.data(), src.data() + src.size(), eig.data());
}

} // namespace msl::eigen_interface

#endif // MSL_EIGEN_INTERFACE_HPP
