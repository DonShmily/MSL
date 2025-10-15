/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \src\matrix\real_matrix_view.cpp
** -----
** File Created: Saturday, 11th October 2025 22:09:45
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Saturday, 11th October 2025 22:49:30
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/


#include "matrix/real_matrix_view.hpp"

namespace msl
{

real_matrix_view real_matrix_view::subview(size_t row_start,
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
    double *sub_data = this->data() + col_start * this->rows_;
    return real_matrix_view(sub_data, this->rows_, sub_cols);
}

real_matrix_view real_matrix_view::column_view(size_t j)
{
    auto col_span = this->column(j);
    return real_matrix_view(col_span.data(), this->rows_, 1);
}

const_real_matrix_view real_matrix_view::column_view(size_t j) const
{
    auto col_span = this->column(j);
    return const_real_matrix_view(col_span.data(), this->rows_, 1);
}

// ---In-place operations ---

real_matrix_view &real_matrix_view::operator+=(const real_matrix_base &other)
{
    assert(this->rows_ == other.rows() && this->cols_ == other.cols());
    for (size_t i = 0; i < this->size(); ++i)
    {
        this->data_[i] += other.data()[i];
    }
    return *this;
}

real_matrix_view &real_matrix_view::operator-=(const real_matrix_base &other)
{
    assert(this->rows_ == other.rows() && this->cols_ == other.cols());
    for (size_t i = 0; i < this->size(); ++i)
    {
        this->data_[i] -= other.data()[i];
    }
    return *this;
}

real_matrix_view &real_matrix_view::operator*=(const double &scalar)
{
    for (size_t i = 0; i < this->size(); ++i)
    {
        this->data_[i] *= scalar;
    }
    return *this;
}

real_matrix_view &real_matrix_view::operator/=(const double &scalar)
{
    for (size_t i = 0; i < this->size(); ++i)
    {
        this->data_[i] /= scalar;
    }
    return *this;
}

// --- Helper functions ---

void swap(real_matrix_view &a, real_matrix_view &b) noexcept { a.swap(b); }

// Create view from owned matrix
[[nodiscard]] real_matrix_view make_view(real_matrix_owned &mat)
{
    return real_matrix_view(mat);
}

[[nodiscard]] const_real_matrix_view make_view(const real_matrix_owned &mat)
{
    return const_real_matrix_view(mat);
}

} // namespace msl