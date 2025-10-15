/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \include\matrix\matrix_view.hpp
** -----
** File Created: Thursday, 9th October 2025 21:37:02
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Friday, 10th October 2025 17:41:37
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_MATRIX_VIEW_HPP
#define MSL_MATRIX_VIEW_HPP

#include <cassert>
#include "matrix_base.hpp"
#include "matrix_owned.hpp"

namespace msl
{

/**
 * @brief Non-owning matrix view
 *
 * Provides a view into existing matrix data without owning it.
 * Lifetime: The viewed data must outlive this view.
 *
 * Copy semantics: Shallow copy (copies the view, not the data)
 *
 * Note: Be careful with lifetime management. This class does not
 * prevent dangling references.
 */
template <typename T>
class const_matrix_view;

template <typename T>
class matrix_view : public matrix_base<T, matrix_view<T>>
{
private:
    using Base = matrix_base<T, matrix_view<T>>;

public:
    // --- Constructors ---

    // Default constructor (empty view)
    matrix_view() : Base() {}

    // Construct from raw pointer
    matrix_view(T *data, size_t rows, size_t cols)
        : Base(std::span<T>(data, rows * cols), rows, cols)
    {}

    // Construct from matrix_owned (most common use case)
    matrix_view(matrix_owned<T> &owner)
        : Base(owner.data_span(), owner.rows(), owner.cols())
    {}

    // Construct from span
    matrix_view(std::span<T> data, size_t rows, size_t cols)
        : Base(data, rows, cols)
    {
        assert(data.size() == rows * cols);
    }

    // --- Copy semantics: Shallow copy (copies the view) ---

    matrix_view(const matrix_view &other) = default;
    matrix_view &operator=(const matrix_view &other) = default;

    // --- Move semantics ---

    matrix_view(matrix_view &&other) noexcept = default;
    matrix_view &operator=(matrix_view &&other) noexcept = default;

    // --- Cannot construct view from const owned (use const_matrix_view) ---
    matrix_view(const matrix_owned<T> &) = delete;

    // --- Destructor ---
    ~matrix_view() = default;

    // --- Subview creation ---

    /**
     * @brief Create a subview (view into a submatrix)
     *
     * Note: The lifetime of the returned view is tied to this view's data
     */
    [[nodiscard]] matrix_view<T>
    subview(size_t row_start, size_t row_end, size_t col_start, size_t col_end)
    {
        assert(row_start < row_end && row_end <= this->rows_);
        assert(col_start < col_end && col_end <= this->cols_);

        // For column-major layout, creating arbitrary subviews requires
        // either strided access or copying. Here we document the limitation.
        // For now, we only support full-column subviews efficiently.
        assert(
            row_start == 0 && row_end == this->rows_
            && "Arbitrary row subviews not supported in column-major layout");

        size_t sub_cols = col_end - col_start;
        T *sub_data = this->data() + col_start * this->rows_;
        return matrix_view<T>(sub_data, this->rows_, sub_cols);
    }

    // --- Column subview (efficient) ---

    [[nodiscard]] matrix_view<T> column_view(size_t j)
    {
        auto col_span = this->column(j);
        return matrix_view<T>(col_span.data(), this->rows_, 1);
    }

    [[nodiscard]] const_matrix_view<T> column_view(size_t j) const
    {
        auto col_span = this->column(j);
        return const_matrix_view<T>(col_span.data(), this->rows_, 1);
    }

    // --- Fill operations ---

    void fill(const T &value) { std::fill(this->begin(), this->end(), value); }

    void fill_zeros() { std::fill(this->begin(), this->end(), T{}); }

    // --- In-place operations ---

    matrix_view &operator+=(const matrix_base<T, matrix_view<T>> &other)
    {
        assert(this->rows_ == other.rows() && this->cols_ == other.cols());
        for (size_t i = 0; i < this->size(); ++i)
        {
            this->data_[i] += other.data()[i];
        }
        return *this;
    }

    matrix_view &operator+=(const matrix_base<T, matrix_owned<T>> &other)
    {
        assert(this->rows_ == other.rows() && this->cols_ == other.cols());
        for (size_t i = 0; i < this->size(); ++i)
        {
            this->data_[i] += other.data()[i];
        }
        return *this;
    }

    matrix_view &operator-=(const matrix_base<T, matrix_view<T>> &other)
    {
        assert(this->rows_ == other.rows() && this->cols_ == other.cols());
        for (size_t i = 0; i < this->size(); ++i)
        {
            this->data_[i] -= other.data()[i];
        }
        return *this;
    }

    matrix_view &operator-=(const matrix_base<T, matrix_owned<T>> &other)
    {
        assert(this->rows_ == other.rows() && this->cols_ == other.cols());
        for (size_t i = 0; i < this->size(); ++i)
        {
            this->data_[i] -= other.data()[i];
        }
        return *this;
    }

    matrix_view &operator*=(const T &scalar)
    {
        for (size_t i = 0; i < this->size(); ++i)
        {
            this->data_[i] *= scalar;
        }
        return *this;
    }

    matrix_view &operator/=(const T &scalar)
    {
        for (size_t i = 0; i < this->size(); ++i)
        {
            this->data_[i] /= scalar;
        }
        return *this;
    }

    // --- Copy data from another matrix ---

    void copy_from(const matrix_base<T, matrix_owned<T>> &src)
    {
        assert(this->rows_ == src.rows() && this->cols_ == src.cols());
        std::copy(src.begin(), src.end(), this->begin());
    }

    void copy_from(const matrix_base<T, matrix_view<T>> &src)
    {
        assert(this->rows_ == src.rows() && this->cols_ == src.cols());
        std::copy(src.begin(), src.end(), this->begin());
    }

    // --- Convert to owned matrix (deep copy) ---

    [[nodiscard]] matrix_owned<T> to_owned() const
    {
        return matrix_owned<T>(this->rows_, this->cols_, this->data_);
    }

    // --- Swap (swaps the views, not the data) ---

    void swap(matrix_view &other) noexcept
    {
        using std::swap;
        swap(this->data_, other.data_);
        swap(this->rows_, other.rows_);
        swap(this->cols_, other.cols_);
    }
};

/**
 * @brief Const matrix view
 *
 * Similar to matrix_view but provides read-only access
 */
template <typename T>
class const_matrix_view
{
private:
    std::span<const T> data_;
    size_t rows_{0};
    size_t cols_{0};

public:
    using value_type = T;
    using const_reference = const T &;
    using const_pointer = const T *;
    using size_type = size_t;

    // --- Constructors ---

    const_matrix_view() = default;

    const_matrix_view(const T *data, size_t rows, size_t cols)
        : data_(data, rows * cols), rows_(rows), cols_(cols)
    {}

    const_matrix_view(const matrix_owned<T> &owner)
        : data_(owner.data_span()), rows_(owner.rows()), cols_(owner.cols())
    {}

    const_matrix_view(const matrix_view<T> &view)
        : data_(view.data_span()), rows_(view.rows()), cols_(view.cols())
    {}

    const_matrix_view(std::span<const T> data, size_t rows, size_t cols)
        : data_(data), rows_(rows), cols_(cols)
    {
        assert(data.size() == rows * cols);
    }

    // --- Copy and move ---

    const_matrix_view(const const_matrix_view &) = default;
    const_matrix_view &operator=(const const_matrix_view &) = default;
    const_matrix_view(const_matrix_view &&) noexcept = default;
    const_matrix_view &operator=(const_matrix_view &&) noexcept = default;

    // --- Size queries ---

    [[nodiscard]] size_t rows() const noexcept { return rows_; }
    [[nodiscard]] size_t cols() const noexcept { return cols_; }
    [[nodiscard]] size_t size() const noexcept { return rows_ * cols_; }
    [[nodiscard]] std::pair<size_t, size_t> shape() const noexcept
    {
        return {rows_, cols_};
    }

    // --- Element access (read-only) ---

    [[nodiscard]] const T &operator()(size_t i, size_t j) const noexcept
    {
        MSL_CHECK_BOUNDS(i < rows_ && j < cols_, "Matrix index out of bounds");
        return data_[j * rows_ + i];
    }

    [[nodiscard]] const T &operator[](size_t idx) const noexcept
    {
        MSL_CHECK_BOUNDS(idx < size(), "Linear index out of bounds");
        return data_[idx];
    }

    // --- Column access ---

    [[nodiscard]] std::span<const T> column(size_t j) const noexcept
    {
        MSL_CHECK_BOUNDS(j < cols_, "Column index out of bounds");
        return {data_.data() + j * rows_, rows_};
    }

    // --- Data access ---

    [[nodiscard]] const T *data() const noexcept { return data_.data(); }
    [[nodiscard]] std::span<const T> data_span() const noexcept
    {
        return data_;
    }

    // --- Iterators ---

    [[nodiscard]] const T *begin() const noexcept { return data_.begin(); }
    [[nodiscard]] const T *end() const noexcept { return data_.end(); }
    [[nodiscard]] const T *cbegin() const noexcept { return data_.begin(); }
    [[nodiscard]] const T *cend() const noexcept { return data_.end(); }

    // --- Convert to owned matrix ---

    [[nodiscard]] matrix_owned<T> to_owned() const
    {
        return matrix_owned<T>(
            this->rows_, this->cols_, std::span<const T>(this->data_));
    }
};

// --- Helper functions ---

template <typename T>
void swap(matrix_view<T> &a, matrix_view<T> &b) noexcept
{
    a.swap(b);
}

// Create view from owned matrix
template <typename T>
[[nodiscard]] matrix_view<T> make_view(matrix_owned<T> &mat)
{
    return matrix_view<T>(mat);
}

template <typename T>
[[nodiscard]] const_matrix_view<T> make_view(const matrix_owned<T> &mat)
{
    return const_matrix_view<T>(mat);
}

} // namespace msl

#endif // MSL_MATRIX_VIEW_HPP