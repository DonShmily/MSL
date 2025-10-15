/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \include\matrix\matrix_owned.hpp
** -----
** File Created: Thursday, 9th October 2025 21:35:42
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Friday, 10th October 2025 17:41:30
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_MATRIX_OWNED_HPP
#define MSL_MATRIX_OWNED_HPP

#include <algorithm>
#include <initializer_list>
#include <vector>
#include "matrix_base.hpp"

namespace msl
{

template <typename T>
class matrix_view; // Forward declaration

/**
 * @brief Matrix with owning semantics
 *
 * Owns its data in a std::vector. Provides deep copy semantics.
 * The span in base class always points to the internal vector.
 */
template <typename T>
class matrix_owned : public matrix_base<T, matrix_owned<T>>
{
private:
    using Base = matrix_base<T, matrix_owned<T>>;
    std::vector<T> storage_; // Actual data storage

    // Update base span when storage changes
    void update_span() { this->data_ = std::span<T>(storage_); }

public:
    // --- Constructors ---

    // Default constructor
    matrix_owned() : Base() {}

    // Size constructor (zero-initialized)
    matrix_owned(size_t rows, size_t cols) : Base(), storage_(rows * cols)
    {
        this->rows_ = rows;
        this->cols_ = cols;
        update_span();
    }

    // Size + value constructor
    matrix_owned(size_t rows, size_t cols, const T &value)
        : Base(), storage_(rows * cols, value)
    {
        this->rows_ = rows;
        this->cols_ = cols;
        update_span();
    }

    // Construct from span (deep copy)
    matrix_owned(size_t rows, size_t cols, std::span<const T> data)
        : Base(), storage_(data.begin(), data.end())
    {
        assert(data.size() == rows * cols);
        this->rows_ = rows;
        this->cols_ = cols;
        update_span();
    }

    // Construct from initializer list (row-major input)
    matrix_owned(size_t rows, size_t cols, std::initializer_list<T> values)
        : Base(), storage_(rows * cols)
    {
        assert(values.size() == rows * cols);
        this->rows_ = rows;
        this->cols_ = cols;

        // Convert row-major to column-major
        auto it = values.begin();
        for (size_t i = 0; i < rows; ++i)
        {
            for (size_t j = 0; j < cols; ++j)
            {
                storage_[j * rows + i] = *it++;
            }
        }
        update_span();
    }

    // --- Copy operations (deep copy) ---

    matrix_owned(const matrix_owned &other) : Base(), storage_(other.storage_)
    {
        this->rows_ = other.rows_;
        this->cols_ = other.cols_;
        update_span();
    }

    matrix_owned &operator=(const matrix_owned &other)
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

    matrix_owned &operator=(const matrix_view<T> &view)
    {
        this->resize(view.rows(), view.cols());
        std::copy(view.begin(), view.end(), this->begin());
        return *this;
    }

    // --- Move operations (no-throw guarantee) ---

    matrix_owned(matrix_owned &&other) noexcept
        : Base(), storage_(std::move(other.storage_))
    {
        this->rows_ = other.rows_;
        this->cols_ = other.cols_;
        update_span();

        other.rows_ = 0;
        other.cols_ = 0;
        other.data_ = {};
    }

    matrix_owned &operator=(matrix_owned &&other) noexcept
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

    // --- Construct from view (deep copy) ---

    explicit matrix_owned(const matrix_view<T> &view)
        : Base(), storage_(view.data(), view.data() + view.size())
    {
        this->rows_ = view.rows();
        this->cols_ = view.cols();
        update_span();
    }

    // --- Destructor ---
    ~matrix_owned() = default;

    // --- Size management ---

    void resize(size_t rows, size_t cols, const T &value = T())
    {
        this->rows_ = rows;
        this->cols_ = cols;
        storage_.resize(rows * cols, value);
        update_span();
    }

    void clear()
    {
        this->rows_ = 0;
        this->cols_ = 0;
        storage_.clear();
        update_span();
    }

    // Reserve capacity (useful for avoiding reallocations)
    void reserve(size_t capacity) { storage_.reserve(capacity); }

    [[nodiscard]] size_t capacity() const noexcept
    {
        return storage_.capacity();
    }

    void shrink_to_fit() { storage_.shrink_to_fit(); }

    // --- Fill operations ---

    void fill(const T &value)
    {
        std::fill(storage_.begin(), storage_.end(), value);
    }

    void fill_zeros() { std::fill(storage_.begin(), storage_.end(), T{}); }

    // --- Row/Column extraction (returns new matrix) ---

    [[nodiscard]] matrix_owned<T> get_row(size_t i) const
    {
        assert(i < this->rows_);
        matrix_owned<T> row(1, this->cols_);
        for (size_t j = 0; j < this->cols_; ++j)
        {
            row(0, j) = (*this)(i, j);
        }
        return row;
    }

    [[nodiscard]] matrix_owned<T> get_column(size_t j) const
    {
        assert(j < this->cols_);
        matrix_owned<T> col(this->rows_, 1);
        const auto &col_span = this->column(j);
        std::copy(col_span.begin(), col_span.end(), col.data());
        return col;
    }

    // --- Submatrix extraction ---

    [[nodiscard]] matrix_owned<T> submatrix(size_t row_start,
                                            size_t row_end,
                                            size_t col_start,
                                            size_t col_end) const
    {
        assert(row_start < row_end && row_end <= this->rows_);
        assert(col_start < col_end && col_end <= this->cols_);

        size_t sub_rows = row_end - row_start;
        size_t sub_cols = col_end - col_start;
        matrix_owned<T> result(sub_rows, sub_cols);

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

    matrix_owned &operator+=(const matrix_base<T, matrix_owned<T>> &other)
    {
        assert(this->rows_ == other.rows() && this->cols_ == other.cols());
        for (size_t i = 0; i < storage_.size(); ++i)
        {
            storage_[i] += other.data()[i];
        }
        return *this;
    }

    matrix_owned &operator-=(const matrix_base<T, matrix_owned<T>> &other)
    {
        assert(this->rows_ == other.rows() && this->cols_ == other.cols());
        for (size_t i = 0; i < storage_.size(); ++i)
        {
            storage_[i] -= other.data()[i];
        }
        return *this;
    }

    matrix_owned &operator*=(const T &scalar)
    {
        for (auto &val : storage_)
        {
            val *= scalar;
        }
        return *this;
    }

    matrix_owned &operator/=(const T &scalar)
    {
        for (auto &val : storage_)
        {
            val /= scalar;
        }
        return *this;
    }

    // --- Swap ---

    void swap(matrix_owned &other) noexcept
    {
        using std::swap;
        swap(storage_, other.storage_);
        swap(this->rows_, other.rows_);
        swap(this->cols_, other.cols_);
        update_span();
        other.update_span();
    }

    // --- Factory methods ---

    [[nodiscard]] static matrix_owned<T> zeros(size_t rows, size_t cols)
    {
        return matrix_owned<T>(rows, cols, T{});
    }

    [[nodiscard]] static matrix_owned<T> ones(size_t rows, size_t cols)
    {
        return matrix_owned<T>(rows, cols, T{1});
    }

    [[nodiscard]] static matrix_owned<T> identity(size_t n)
    {
        matrix_owned<T> result(n, n, T{});
        for (size_t i = 0; i < n; ++i)
        {
            result(i, i) = T{1};
        }
        return result;
    }

    // Create diagonal matrix from vector
    [[nodiscard]] static matrix_owned<T> diagonal(std::span<const T> diag)
    {
        size_t n = diag.size();
        matrix_owned<T> result(n, n, T{});
        for (size_t i = 0; i < n; ++i)
        {
            result(i, i) = diag[i];
        }
        return result;
    }
};

// --- Free functions ---

template <typename T>
void swap(matrix_owned<T> &a, matrix_owned<T> &b) noexcept
{
    a.swap(b);
}

// Binary operators (return new matrix)
template <typename T>
[[nodiscard]] matrix_owned<T> operator+(const matrix_owned<T> &a,
                                        const matrix_owned<T> &b)
{
    assert(a.rows() == b.rows() && a.cols() == b.cols());
    matrix_owned<T> result(a);
    result += b;
    return result;
}

template <typename T>
[[nodiscard]] matrix_owned<T> operator-(const matrix_owned<T> &a,
                                        const matrix_owned<T> &b)
{
    assert(a.rows() == b.rows() && a.cols() == b.cols());
    matrix_owned<T> result(a);
    result -= b;
    return result;
}

template <typename T>
[[nodiscard]] matrix_owned<T> operator*(const matrix_owned<T> &mat,
                                        const T &scalar)
{
    matrix_owned<T> result(mat);
    result *= scalar;
    return result;
}

template <typename T>
[[nodiscard]] matrix_owned<T> operator*(const T &scalar,
                                        const matrix_owned<T> &mat)
{
    return mat * scalar;
}

template <typename T>
[[nodiscard]] matrix_owned<T> operator/(const matrix_owned<T> &mat,
                                        const T &scalar)
{
    matrix_owned<T> result(mat);
    result /= scalar;
    return result;
}

} // namespace msl

#endif // MSL_MATRIX_OWNED_HPP