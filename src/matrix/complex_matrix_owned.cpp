/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \src\matrix\complex_matrix_owned.cpp
** -----
** File Created: Saturday, 11th October 2025 22:09:45
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Saturday, 11th October 2025 22:43:08
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#include <algorithm>
#include <cassert>

#include "matrix/complex_matrix_owned.hpp"
#include "matrix/complex_matrix_view.hpp"

namespace msl
{
// --- Copy/Move Operations ---

complex_matrix_owned::complex_matrix_owned(const complex_matrix_owned &other)
    : Base(), storage_(other.storage_)
{
    this->rows_ = other.rows_;
    this->cols_ = other.cols_;
    update_span();
}

complex_matrix_owned::complex_matrix_owned(const complex_matrix_view &view)
    : Base()
{
    this->rows_ = view.rows();
    this->cols_ = view.cols();
    storage_.resize(view.size());
    // Manual copy required as view might not be contiguous column-major
    for (size_t j = 0; j < view.cols(); ++j)
    {
        for (size_t i = 0; i < view.rows(); ++i)
        {
            (*this)(i, j) = view(i, j);
        }
    }
    update_span();
}

complex_matrix_owned &
complex_matrix_owned::operator=(const complex_matrix_owned &other)
{
    if (this != &other)
    {
        storage_ = other.storage_;
        this->rows_ = other.rows_;
        this->cols_ = other.cols_;
        update_span();
    }
    return *this;
}

complex_matrix_owned &
complex_matrix_owned::operator=(const complex_matrix_view &view)
{
    this->resize(view.rows(), view.cols());
    // Note: std::copy cannot be used directly if view is not contiguous
    // column-major. A manual loop is safer.
    for (size_t j = 0; j < view.cols(); ++j)
    {
        for (size_t i = 0; i < view.rows(); ++i)
        {
            (*this)(i, j) = view(i, j);
        }
    }
    update_span();
    return *this;
}

complex_matrix_owned::complex_matrix_owned(
    complex_matrix_owned &&other) noexcept
    : Base(), storage_(std::move(other.storage_))
{
    this->rows_ = other.rows_;
    this->cols_ = other.cols_;
    update_span();

    other.rows_ = 0;
    other.cols_ = 0;
    other.data_ = {};
}

complex_matrix_owned &
complex_matrix_owned::operator=(complex_matrix_owned &&other) noexcept
{
    if (this != &other)
    {
        storage_ = std::move(other.storage_);
        this->rows_ = other.rows_;
        this->cols_ = other.cols_;
        update_span();

        other.rows_ = 0;
        other.cols_ = 0;
        other.data_ = {};
    }
    return *this;
}

// --- Row/Column extraction ---

complex_matrix_owned complex_matrix_owned::get_row(size_t i) const
{
    assert(i < this->rows_);
    complex_matrix_owned row(1, this->cols_);
    for (size_t j = 0; j < this->cols_; ++j)
    {
        row(0, j) = (*this)(i, j);
    }
    return row;
}

complex_matrix_owned complex_matrix_owned::get_column(size_t j) const
{
    assert(j < this->cols_);
    complex_matrix_owned col(this->rows_, 1);
    const auto &col_span = this->column(j);
    std::copy(col_span.begin(), col_span.end(), col.storage_.begin());
    return col;
}

// --- Submatrix extraction ---

complex_matrix_owned complex_matrix_owned::submatrix(size_t row_start,
                                                     size_t row_end,
                                                     size_t col_start,
                                                     size_t col_end) const
{
    assert(row_start < row_end && row_end <= this->rows_);
    assert(col_start < col_end && col_end <= this->cols_);

    size_t sub_rows = row_end - row_start;
    size_t sub_cols = col_end - col_start;
    complex_matrix_owned result(sub_rows, sub_cols);

    for (size_t j = 0; j < sub_cols; ++j)
    {
        for (size_t i = 0; i < sub_rows; ++i)
        {
            result(i, j) = (*this)(row_start + i, col_start + j);
        }
    }
    return result;
}

// --- In-place operations ---

complex_matrix_owned &
complex_matrix_owned::operator+=(const complex_matrix_base &other)
{
    assert(this->rows_ == other.rows() && this->cols_ == other.cols());
    for (size_t i = 0; i < storage_.size(); ++i)
    {
        storage_[i] += other.data()[i];
    }
    return *this;
}

complex_matrix_owned &
complex_matrix_owned::operator-=(const complex_matrix_base &other)
{
    assert(this->rows_ == other.rows() && this->cols_ == other.cols());
    for (size_t i = 0; i < storage_.size(); ++i)
    {
        storage_[i] -= other.data()[i];
    }
    return *this;
}

complex_matrix_owned &
complex_matrix_owned::operator*=(const std::complex<double> &scalar)
{
    for (auto &val : storage_)
    {
        val *= scalar;
    }
    return *this;
}
complex_matrix_owned &
complex_matrix_owned::operator/=(const std::complex<double> &scalar)
{
    for (auto &val : storage_)
    {
        val /= scalar;
    }
    return *this;
}

// --- Apply function ---
void complex_matrix_owned::apply(
    std::function<std::complex<double>(std::complex<double>)> func)
{
    for (auto &val : storage_)
    {
        val = func(val);
    }
}

// --- Swap ---

void complex_matrix_owned::swap(complex_matrix_owned &other) noexcept
{
    using std::swap;
    swap(storage_, other.storage_);
    swap(this->rows_, other.rows_);
    swap(this->cols_, other.cols_);
    update_span();
    other.update_span();
}

// --- Factory methods ---
complex_matrix_owned
complex_matrix_owned::diagonal(std::span<const std::complex<double>> diag)
{
    size_t n = diag.size();
    complex_matrix_owned result(n, n, 0.0);
    for (size_t i = 0; i < n; ++i)
    {
        result(i, i) = diag[i];
    }
    return result;
}

// --- Free functions ---
void swap(complex_matrix_owned &a, complex_matrix_owned &b) noexcept
{
    a.swap(b);
}

// Binary operators (return new matrix)
[[nodiscard]] complex_matrix_owned operator+(const complex_matrix_owned &a,
                                             const complex_matrix_owned &b)
{
    assert(a.rows() == b.rows() && a.cols() == b.cols());
    complex_matrix_owned result(a);
    result += b;
    return result;
}

[[nodiscard]] complex_matrix_owned operator-(const complex_matrix_owned &a,
                                             const complex_matrix_owned &b)
{
    assert(a.rows() == b.rows() && a.cols() == b.cols());
    complex_matrix_owned result(a);
    result -= b;
    return result;
}

[[nodiscard]] complex_matrix_owned operator*(const complex_matrix_owned &mat,
                                             const double &scalar)
{
    complex_matrix_owned result(mat);
    result *= scalar;
    return result;
}

[[nodiscard]] complex_matrix_owned operator*(const double &scalar,
                                             const complex_matrix_owned &mat)
{
    return mat * scalar;
}

[[nodiscard]] complex_matrix_owned operator/(const complex_matrix_owned &mat,
                                             const double &scalar)
{
    complex_matrix_owned result(mat);
    result /= scalar;
    return result;
}

} // namespace msl