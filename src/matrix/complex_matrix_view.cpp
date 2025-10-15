/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \src\matrix\complex_matrix_view.cpp
** -----
** File Created: Saturday, 11th October 2025 22:09:45
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Saturday, 11th October 2025 22:49:30
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#include "matrix/complex_matrix_view.hpp"

namespace msl
{

complex_matrix_view complex_matrix_view::subview(size_t row_start,
                                                 size_t row_end,
                                                 size_t col_start,
                                                 size_t col_end)
{
    assert(row_start < row_end && row_end <= this->rows_);
    assert(col_start < col_end && col_end <= this->cols_);

    // For column-major layout, creating arbitrary subviews requires
    // either strided access or copying. Here we document the limitation.
    // For now, we only support full-column subviews efficiently.
    assert(row_start == 0 && row_end == this->rows_
           && "Arbitrary row subviews not supported in column-major layout");

    size_t sub_cols = col_end - col_start;
    std::complex<double> *sub_data = this->data() + col_start * this->rows_;
    return complex_matrix_view(sub_data, this->rows_, sub_cols);
}

complex_matrix_view complex_matrix_view::column_view(size_t j)
{
    auto col_span = this->column(j);
    return complex_matrix_view(col_span.data(), this->rows_, 1);
}

const_complex_matrix_view complex_matrix_view::column_view(size_t j) const
{
    auto col_span = this->column(j);
    return const_complex_matrix_view(col_span.data(), this->rows_, 1);
}

// ---In-place operations ---

complex_matrix_view &
complex_matrix_view::operator+=(const complex_matrix_base &other)
{
    assert(this->rows_ == other.rows() && this->cols_ == other.cols());
    for (size_t i = 0; i < this->size(); ++i)
    {
        this->data_[i] += other.data()[i];
    }
    return *this;
}

complex_matrix_view &
complex_matrix_view::operator-=(const complex_matrix_base &other)
{
    assert(this->rows_ == other.rows() && this->cols_ == other.cols());
    for (size_t i = 0; i < this->size(); ++i)
    {
        this->data_[i] -= other.data()[i];
    }
    return *this;
}

complex_matrix_view &
complex_matrix_view::operator*=(const std::complex<double> &scalar)
{
    for (size_t i = 0; i < this->size(); ++i)
    {
        this->data_[i] *= scalar;
    }
    return *this;
}

complex_matrix_view &
complex_matrix_view::operator/=(const std::complex<double> &scalar)
{
    for (size_t i = 0; i < this->size(); ++i)
    {
        this->data_[i] /= scalar;
    }
    return *this;
}

// --- Helper functions ---

void swap(complex_matrix_view &a, complex_matrix_view &b) noexcept
{
    a.swap(b);
}

// Create view from owned matrix
[[nodiscard]] complex_matrix_view make_view(complex_matrix_owned &mat)
{
    return complex_matrix_view(mat);
}

[[nodiscard]] const_complex_matrix_view
make_view(const complex_matrix_owned &mat)
{
    return const_complex_matrix_view(mat);
}

} // namespace msl