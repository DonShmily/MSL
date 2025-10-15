/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \include\matrix\matrix_real.hpp
** -----
** File Created: Friday, 10th October 2025 14:36:41
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Friday, 10th October 2025 17:41:47
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/


#ifndef MSL_MATRIX_REAL_HPP
#define MSL_MATRIX_REAL_HPP

#include <algorithm>
#include <cmath>
#include <numeric>
#include "matrix_owned.hpp"

namespace msl
{

/**
 * @brief Real-valued matrix (double precision)
 *
 * Specialization of matrix_owned for double with additional
 * real-matrix specific operations.
 */
class matrix_real : public matrix_owned<double>
{
public:
    using value_type = double;
    using Base = matrix_owned<double>;

    // Inherit all constructors from base
    using Base::matrix_owned;

    // --- Construct from Base (matrix_owned) ---
    matrix_real(const matrix_owned<double> &other) : Base(other) {}
    matrix_real(matrix_owned<double> &&other) noexcept : Base(std::move(other))
    {}

    // --- Real matrix specific operations ---

    /**
     * @brief Compute matrix trace (sum of diagonal elements)
     */
    [[nodiscard]] double trace() const
    {
        assert(this->rows_ == this->cols_);
        double sum = 0.0;
        for (size_t i = 0; i < this->rows_; ++i)
        {
            sum += (*this)(i, i);
        }
        return sum;
    }

    /**
     * @brief Compute Frobenius norm (sqrt of sum of squares)
     */
    [[nodiscard]] double norm_frobenius() const
    {
        double sum = 0.0;
        for (const auto &val : *this)
        {
            sum += val * val;
        }
        return std::sqrt(sum);
    }

    /**
     * @brief Compute infinity norm (max absolute row sum)
     */
    [[nodiscard]] double norm_inf() const
    {
        double max_sum = 0.0;
        for (size_t i = 0; i < this->rows_; ++i)
        {
            double row_sum = 0.0;
            for (size_t j = 0; j < this->cols_; ++j)
            {
                row_sum += std::abs((*this)(i, j));
            }
            max_sum = std::max(max_sum, row_sum);
        }
        return max_sum;
    }

    /**
     * @brief Compute 1-norm (max absolute column sum)
     */
    [[nodiscard]] double norm_1() const
    {
        double max_sum = 0.0;
        for (size_t j = 0; j < this->cols_; ++j)
        {
            double col_sum = 0.0;
            auto col = this->column(j);
            for (const auto &val : col)
            {
                col_sum += std::abs(val);
            }
            max_sum = std::max(max_sum, col_sum);
        }
        return max_sum;
    }

    /**
     * @brief Transpose matrix
     */
    [[nodiscard]] matrix_real transpose() const
    {
        matrix_real result(this->cols_, this->rows_);
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
    [[nodiscard]] std::vector<double> diagonal() const
    {
        size_t n = std::min(this->rows_, this->cols_);
        std::vector<double> diag(n);
        for (size_t i = 0; i < n; ++i)
        {
            diag[i] = (*this)(i, i);
        }
        return diag;
    }

    /**
     * @brief Extract lower triangular part
     */
    [[nodiscard]] matrix_real
    lower_triangular(bool include_diagonal = true) const
    {
        assert(this->rows_ == this->cols_);
        matrix_real result(this->rows_, this->cols_, 0.0);

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
    [[nodiscard]] matrix_real
    upper_triangular(bool include_diagonal = true) const
    {
        assert(this->rows_ == this->cols_);
        matrix_real result(this->rows_, this->cols_, 0.0);

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

    /**
     * @brief Element-wise absolute value
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
     * @brief Apply function element-wise
     */
    template <typename Func>
    [[nodiscard]] matrix_real apply(Func &&func) const
    {
        matrix_real result(this->rows_, this->cols_);
        for (size_t i = 0; i < this->size(); ++i)
        {
            result[i] = func((*this)[i]);
        }
        return result;
    }

    /**
     * @brief Sum of all elements
     */
    [[nodiscard]] double sum() const
    {
        return std::accumulate(this->begin(), this->end(), 0.0);
    }

    /**
     * @brief Mean of all elements
     */
    [[nodiscard]] double mean() const
    {
        return sum() / static_cast<double>(this->size());
    }

    /**
     * @brief Minimum element
     */
    [[nodiscard]] double min() const
    {
        return *std::min_element(this->begin(), this->end());
    }

    /**
     * @brief Maximum element
     */
    [[nodiscard]] double max() const
    {
        return *std::max_element(this->begin(), this->end());
    }

    /**
     * @brief Row sums
     */
    [[nodiscard]] std::vector<double> row_sums() const
    {
        std::vector<double> sums(this->rows_, 0.0);
        for (size_t i = 0; i < this->rows_; ++i)
        {
            for (size_t j = 0; j < this->cols_; ++j)
            {
                sums[i] += (*this)(i, j);
            }
        }
        return sums;
    }

    /**
     * @brief Column sums
     */
    [[nodiscard]] std::vector<double> column_sums() const
    {
        std::vector<double> sums(this->cols_, 0.0);
        for (size_t j = 0; j < this->cols_; ++j)
        {
            auto col = this->column(j);
            sums[j] = std::accumulate(col.begin(), col.end(), 0.0);
        }
        return sums;
    }
};

// Type alias for convenience
using matrixd = matrix_real;

// --- Matrix multiplication ---

/**
 * @brief Matrix-matrix multiplication (C = A * B)
 *
 * Note: This is a naive implementation. For production use,
 * integrate with BLAS (dgemm) for better performance.
 */
[[nodiscard]] inline matrix_real operator*(const matrix_real &A,
                                           const matrix_real &B)
{
    assert(A.cols() == B.rows());

    matrix_real C(A.rows(), B.cols(), 0.0);

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
 * @brief Matrix-vector multiplication (y = A * x)
 */
[[nodiscard]] inline std::vector<double> operator*(const matrix_real &A,
                                                   const std::vector<double> &x)
{
    assert(A.cols() == x.size());

    std::vector<double> y(A.rows(), 0.0);

    for (size_t j = 0; j < A.cols(); ++j)
    {
        for (size_t i = 0; i < A.rows(); ++i)
        {
            y[i] += A(i, j) * x[j];
        }
    }

    return y;
}

} // namespace msl

#endif // MSL_MATRIX_REAL_HPP