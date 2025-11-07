/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025
**  File: martrix_decompose.hpp
**  Description: Matrix Decomposition Functions
*/

#ifndef MSL_MATRIX_DECOMPOSE_HPP
#define MSL_MATRIX_DECOMPOSE_HPP

#include <array>
#include <vector>

#include "complex_matrix_base.hpp"
#include "complex_matrix_owned.hpp"
#include "eigen_interface.hpp"
#include "real_matrix_base.hpp"
#include "real_matrix_owned.hpp"

namespace msl::matrix
{
// ============================================================================
// 1. Matrix SVD Functions
// ============================================================================
// Decompose real matrix A into U, S, V such that A = U * S * V^T
// U and V are orthogonal matrices, S is diagonal matrix of singular values
inline std::array<real_matrix_owned, 3> svd(const real_matrix_base &A)
{
    // Placeholder implementation (to be replaced with actual SVD algorithm)
    size_t m = A.rows();
    size_t n = A.cols();
    std::array<real_matrix_owned, 3> result;
    auto &U = result[0] = real_matrix_owned(m, m); // Orthogonal matrix
    auto &S = result[1] = real_matrix_owned(m, n); // Diagonal matrix
    auto &V = result[2] = real_matrix_owned(n, n); // Orthogonal matrix

    auto eig_A = eigen_interface::as_eigen(A);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        eig_A, Eigen::ComputeFullU | Eigen::ComputeFullV);

    auto eig_U = svd.matrixU();
    auto eig_V = svd.matrixV();
    auto eig_S = svd.singularValues();

    // Copy
    std::copy(eig_U.data(), eig_U.data() + eig_U.size(), U.data());
    std::copy(eig_V.data(), eig_V.data() + eig_V.size(), V.data());
    std::fill(S.data(), S.data() + S.size(), 0.0);
    for (size_t i = 0; i < std::min(m, n); ++i)
    {
        S(i, i) = eig_S(i);
    }

    return result;
}

inline std::array<complex_matrix_owned, 3> svd(const complex_matrix_base &A)
{
    // Placeholder implementation (to be replaced with actual SVD algorithm)
    size_t m = A.rows();
    size_t n = A.cols();
    std::array<complex_matrix_owned, 3> result;
    auto &U = result[0] = complex_matrix_owned(m, m); // Orthogonal matrix
    auto &S = result[1] = complex_matrix_owned(m, n); // Diagonal matrix
    auto &V = result[2] = complex_matrix_owned(n, n); // Orthogonal matrix

    auto eig_A = eigen_interface::as_eigen(A);
    Eigen::JacobiSVD<
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>>
        svd(eig_A, Eigen::ComputeFullU | Eigen::ComputeFullV);

    auto eig_U = svd.matrixU();
    auto eig_V = svd.matrixV();
    auto eig_S = svd.singularValues();

    // Copy
    std::copy(eig_U.data(), eig_U.data() + eig_U.size(), U.data());
    std::copy(eig_V.data(), eig_V.data() + eig_V.size(), V.data());
    std::fill(S.data(), S.data() + S.size(), 0.0);
    for (size_t i = 0; i < std::min(m, n); ++i)
    {
        S(i, i) = eig_S(i);
    }

    return result;
}

} // namespace msl::matrix


#endif // MSL_MATRIX_DECOMPOSE_HPP