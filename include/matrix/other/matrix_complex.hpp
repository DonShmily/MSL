/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \include\matrix\matrix_complex.hpp
** -----
** File Created: Friday, 10th October 2025 14:36:35
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Friday, 10th October 2025 17:41:54
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/


#ifndef MSL_MATRIX_COMPLEX_HPP
#define MSL_MATRIX_COMPLEX_HPP

#include <algorithm>
#include <cmath>
#include <complex>
#include <numeric>
#include "matrix_owned.hpp"
#include "matrix_real.hpp"


namespace msl
{

/**
 * @brief Complex-valued matrix
 *
 * Specialization of matrix_owned for std::complex<double> with additional
 * complex-matrix specific operations.
 */
class matrix_complex : public matrix_owned<std::complex<double>>
{
public:
    using value_type = std::complex<double>;
    using Base = matrix_owned<std::complex<double>>;

    // Inherit all constructors from base
    using Base::matrix_owned;

    // --- Construct from real and imaginary parts ---

    matrix_complex(const matrix_real &real_part, const matrix_real &imag_part)
        : Base(real_part.rows(), real_part.cols())
    {
        assert(real_part.rows() == imag_part.rows());
        assert(real_part.cols() == imag_part.cols());

        for (size_t i = 0; i < this->size(); ++i)
        {
            (*this)[i] = std::complex<double>(real_part[i], imag_part[i]);
        }
    }

    // Construct from real matrix (imaginary part = 0)
    explicit matrix_complex(const matrix_real &real_part)
        : Base(real_part.rows(), real_part.cols())
    {
        for (size_t i = 0; i < this->size(); ++i)
        {
            (*this)[i] = std::complex<double>(real_part[i], 0.0);
        }
    }

    // --- Complex matrix specific operations ---

    /**
     * @brief Complex conjugate of the matrix
     */
    [[nodiscard]] matrix_complex conjugate() const
    {
        matrix_complex result(this->rows_, this->cols_);
        for (size_t i = 0; i < this->size(); ++i)
        {
            result[i] = std::conj((*this)[i]);
        }
        return result;
    }

    /**
     * @brief Conjugate transpose (Hermitian transpose)
     */
    [[nodiscard]] matrix_complex conjugate_transpose() const
    {
        matrix_complex result(this->cols_, this->rows_);
        for (size_t i = 0; i < this->rows_; ++i)
        {
            for (size_t j = 0; j < this->cols_; ++j)
            {
                result(j, i) = std::conj((*this)(i, j));
            }
        }
        return result;
    }

    /**
     * @brief Alias for conjugate_transpose
     */
    [[nodiscard]] matrix_complex adjoint() const
    {
        return conjugate_transpose();
    }

    /**
     * @brief Extract real part as real matrix
     */
    [[nodiscard]] matrix_real real() const
    {
        matrix_real result(this->rows_, this->cols_);
        for (size_t i = 0; i < this->size(); ++i)
        {
            result[i] = (*this)[i].real();
        }
        return result;
    }

    /**
     * @brief Extract imaginary part as real matrix
     */
    [[nodiscard]] matrix_real imag() const
    {
        matrix_real result(this->rows_, this->cols_);
        for (size_t i = 0; i < this->size(); ++i)
        {
            result[i] = (*this)[i].imag();
        }
        return result;
    }

    /**
     * @brief Element-wise absolute value (magnitude)
     */
    [[nodiscard]] matrix_real abs() const
    {
        matrix_real result(this->rows_, this->cols_);
        for (size_t i = 0; i < this->size(); ++i)
        {
            result[i] = std::abs((*this)[i]);
        }
        return result;
    }

    /**
     * @brief Element-wise phase (argument)
     */
    [[nodiscard]] matrix_real arg() const
    {
        matrix_real result(this->rows_, this->cols_);
        for (size_t i = 0; i < this->size(); ++i)
        {
            result[i] = std::arg((*this)[i]);
        }
        return result;
    }

    /**
     * @brief Compute matrix trace (sum of diagonal elements)
     */
    [[nodiscard]] std::complex<double> trace() const
    {
        assert(this->rows_ == this->cols_);
        std::complex<double> sum = 0.0;
        for (size_t i = 0; i < this->rows_; ++i)
        {
            sum += (*this)(i, i);
        }
        return sum;
    }

    /**
     * @brief Compute Frobenius norm
     */
    [[nodiscard]] double norm_frobenius() const
    {
        double sum = 0.0;
        for (const auto &val : *this)
        {
            sum += std::norm(val); // |z|^2
        }
        return std::sqrt(sum);
    }

    /**
     * @brief Transpose (without conjugation)
     */
    [[nodiscard]] matrix_complex transpose() const
    {
        matrix_complex result(this->cols_, this->rows_);
        for (size_t i = 0; i < this->rows_; ++i)
        {
            for (size_t j = 0; j < this->cols_; ++j)
            {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
    }

    /**
     * @brief Transpose in place (only for square matrices)
     */
    void transpose_in_place()
    {
        assert(this->rows_ == this->cols_);
        for (size_t i = 0; i < this->rows_; ++i)
        {
            for (size_t j = i + 1; j < this->cols_; ++j)
            {
                std::swap((*this)(i, j), (*this)(j, i));
            }
        }
    }

    /**
     * @brief Get diagonal elements as a vector
     */
    [[nodiscard]] std::vector<std::complex<double>> diagonal() const
    {
        size_t n = std::min(this->rows_, this->cols_);
        std::vector<std::complex<double>> diag(n);
        for (size_t i = 0; i < n; ++i)
        {
            diag[i] = (*this)(i, i);
        }
        return diag;
    }

    /**
     * @brief Apply function element-wise
     */
    template <typename Func>
    [[nodiscard]] matrix_complex apply(Func &&func) const
    {
        matrix_complex result(this->rows_, this->cols_);
        for (size_t i = 0; i < this->size(); ++i)
        {
            result[i] = func((*this)[i]);
        }
        return result;
    }

    /**
     * @brief Sum of all elements
     */
    [[nodiscard]] std::complex<double> sum() const
    {
        return std::accumulate(
            this->begin(), this->end(), std::complex<double>(0.0, 0.0));
    }

    /**
     * @brief Mean of all elements
     */
    [[nodiscard]] std::complex<double> mean() const
    {
        return sum() / static_cast<double>(this->size());
    }

    /**
     * @brief Extract lower triangular part
     */
    [[nodiscard]] matrix_complex
    lower_triangular(bool include_diagonal = true) const
    {
        assert(this->rows_ == this->cols_);
        matrix_complex result(this->rows_, this->cols_, 0.0);

        for (size_t i = 0; i < this->rows_; ++i)
        {
            size_t j_max = include_diagonal ? i + 1 : i;
            for (size_t j = 0; j < j_max; ++j)
            {
                result(i, j) = (*this)(i, j);
            }
        }
        return result;
    }

    /**
     * @brief Extract upper triangular part
     */
    [[nodiscard]] matrix_complex
    upper_triangular(bool include_diagonal = true) const
    {
        assert(this->rows_ == this->cols_);
        matrix_complex result(this->rows_, this->cols_, 0.0);

        for (size_t i = 0; i < this->rows_; ++i)
        {
            size_t j_start = include_diagonal ? i : i + 1;
            for (size_t j = j_start; j < this->cols_; ++j)
            {
                result(i, j) = (*this)(i, j);
            }
        }
        return result;
    }
};

// Type alias for convenience
using matrixc = matrix_complex;

// --- Matrix multiplication ---

/**
 * @brief Complex matrix-matrix multiplication (C = A * B)
 *
 * Note: This is a naive implementation. For production use,
 * integrate with BLAS (zgemm) for better performance.
 */
[[nodiscard]] inline matrix_complex operator*(const matrix_complex &A,
                                              const matrix_complex &B)
{
    assert(A.cols() == B.rows());

    matrix_complex C(A.rows(), B.cols(), 0.0);

    // Column-major friendly order
    for (size_t j = 0; j < B.cols(); ++j)
    {
        for (size_t k = 0; k < A.cols(); ++k)
        {
            for (size_t i = 0; i < A.rows(); ++i)
            {
                C(i, j) += A(i, k) * B(k, j);
            }
        }
    }

    return C;
}

/**
 * @brief Complex matrix-vector multiplication (y = A * x)
 */
[[nodiscard]] inline std::vector<std::complex<double>>
operator*(const matrix_complex &A, const std::vector<std::complex<double>> &x)
{
    assert(A.cols() == x.size());

    std::vector<std::complex<double>> y(A.rows(), 0.0);

    for (size_t j = 0; j < A.cols(); ++j)
    {
        for (size_t i = 0; i < A.rows(); ++i)
        {
            y[i] += A(i, j) * x[j];
        }
    }

    return y;
}

/**
 * @brief Real matrix * complex matrix
 */
[[nodiscard]] inline matrix_complex operator*(const matrix_real &A,
                                              const matrix_complex &B)
{
    return matrix_complex(A) * B;
}

/**
 * @brief Complex matrix * real matrix
 */
[[nodiscard]] inline matrix_complex operator*(const matrix_complex &A,
                                              const matrix_real &B)
{
    return A * matrix_complex(B);
}

} // namespace msl

#endif // MSL_MATRIX_COMPLEX_HPP