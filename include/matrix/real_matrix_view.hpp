/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \include\matrix\real_matrix_view.hpp
** -----
** File Created: Saturday, 11th October 2025 21:20:53
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Saturday, 11th October 2025 22:49:24
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/


#ifndef MSL_REAL_MATRIX_VIEW_HPP
#define MSL_REAL_MATRIX_VIEW_HPP

#include <cassert>
#include "real_matrix_base.hpp"
#include "real_matrix_owned.hpp"

namespace msl
{
class const_real_matrix_view; // Forward declaration

class real_matrix_view : public real_matrix_base
{
private:
    using Base = real_matrix_base;

public:
    // --- Constructors ---
    // Default constructor (empty view)
    real_matrix_view() : Base() {}

    // Construct from raw pointer
    real_matrix_view(double *data, size_t rows, size_t cols)
        : Base(rows, cols, std::span<double>(data, rows * cols))
    {}

    // Construct from real_matrix_owned (most common use case)
    real_matrix_view(real_matrix_owned &owner)
        : Base(owner.rows(), owner.cols(), owner.data_span())
    {}

    // Construct from span
    real_matrix_view(std::span<double> data, size_t rows, size_t cols)
        : Base(rows, cols, data)
    {
        assert(data.size() == rows * cols);
    }

    // --- Copy semantics: Shallow copy (copies the view) ---
    real_matrix_view(const real_matrix_view &other) = default;
    real_matrix_view &operator=(const real_matrix_view &other) = default;

    // --- Move semantics ---
    real_matrix_view(real_matrix_view &&other) noexcept = default;
    real_matrix_view &operator=(real_matrix_view &&other) noexcept = default;

    // --- Cannot construct view from const owned (use const_matrix_view) ---
    real_matrix_view(const real_matrix_owned &) = delete;

    // --- Destructor ---
    ~real_matrix_view() = default;

    // --- Subview creation ---
    [[nodiscard]] real_matrix_view
    subview(size_t row_start, size_t row_end, size_t col_start, size_t col_end);

    // --- Column subview (efficient) ---
    [[nodiscard]] real_matrix_view column_view(size_t j);
    [[nodiscard]] const_real_matrix_view column_view(size_t j) const;

    // --- Fill operations ---
    void fill(const double &value)
    {
        std::fill(this->begin(), this->end(), value);
    }

    void fill_zeros() { std::fill(this->begin(), this->end(), 0.0); }

    // --- In-place operations ---
    real_matrix_view &operator+=(const real_matrix_base &other);
    real_matrix_view &operator-=(const real_matrix_base &other);
    real_matrix_view &operator*=(const double &scalar);
    real_matrix_view &operator/=(const double &scalar);

    // --- Copy data from another matrix ---
    void copy_from(const real_matrix_base &src)
    {
        assert(this->rows_ == src.rows() && this->cols_ == src.cols());
        std::copy(src.begin(), src.end(), this->begin());
    }

    // --- Convert to owned matrix (deep copy) ---
    [[nodiscard]] real_matrix_owned to_owned() const
    {
        return real_matrix_owned(this->rows_, this->cols_, this->data_);
    }

    // --- Swap (swaps the views, not the data) ---
    void swap(real_matrix_view &other) noexcept
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
class const_real_matrix_view
{
private:
    std::span<const double> data_;
    size_t rows_{0};
    size_t cols_{0};

public:
    // --- Constructors ---
    const_real_matrix_view() = default;

    const_real_matrix_view(const double *data, size_t rows, size_t cols)
        : data_(data, rows * cols), rows_(rows), cols_(cols)
    {}

    const_real_matrix_view(const real_matrix_owned &owner)
        : data_(owner.data_span()), rows_(owner.rows()), cols_(owner.cols())
    {}

    const_real_matrix_view(const real_matrix_view &view)
        : data_(view.data_span()), rows_(view.rows()), cols_(view.cols())
    {}

    const_real_matrix_view(std::span<const double> data,
                           size_t rows,
                           size_t cols)
        : data_(data), rows_(rows), cols_(cols)
    {
        assert(data.size() == rows * cols);
    }

    // --- Copy and move ---
    const_real_matrix_view(const const_real_matrix_view &) = default;
    const_real_matrix_view &operator=(const const_real_matrix_view &) = default;
    const_real_matrix_view(const_real_matrix_view &&) noexcept = default;
    const_real_matrix_view &
    operator=(const_real_matrix_view &&) noexcept = default;

    // --- Size queries ---
    [[nodiscard]] size_t rows() const noexcept { return rows_; }
    [[nodiscard]] size_t cols() const noexcept { return cols_; }
    [[nodiscard]] size_t size() const noexcept { return rows_ * cols_; }
    [[nodiscard]] std::pair<size_t, size_t> shape() const noexcept
    {
        return {rows_, cols_};
    }

    // --- Element access (read-only) ---
    [[nodiscard]] const double &operator()(size_t i, size_t j) const noexcept
    {
        return data_[j * rows_ + i];
    }

    [[nodiscard]] const double &operator[](size_t idx) const noexcept
    {
        return data_[idx];
    }

    // --- Column access ---
    [[nodiscard]] std::span<const double> column(size_t j) const noexcept
    {
        return {data_.data() + j * rows_, rows_};
    }

    // --- Data access ---
    [[nodiscard]] const double *data() const noexcept { return data_.data(); }
    [[nodiscard]] std::span<const double> data_span() const noexcept
    {
        return data_;
    }

    // --- Iterators ---
    [[nodiscard]] auto begin() const noexcept { return data_.begin(); }
    [[nodiscard]] auto end() const noexcept { return data_.end(); }
    [[nodiscard]] auto cbegin() const noexcept { return data_.begin(); }
    [[nodiscard]] auto cend() const noexcept { return data_.end(); }

    // --- Convert to owned matrix ---
    [[nodiscard]] real_matrix_owned to_owned() const
    {
        return real_matrix_owned(
            this->rows_, this->cols_, std::span<const double>(this->data_));
    }
};
typedef real_matrix_view matrixd_view;
typedef const_real_matrix_view const_matrixd_view;

// --- Helper functions ---
void swap(real_matrix_view &a, real_matrix_view &b) noexcept;
// Create view from owned matrix
[[nodiscard]] real_matrix_view make_view(real_matrix_owned &mat);
[[nodiscard]] const_real_matrix_view make_view(const real_matrix_owned &mat);

} // namespace msl

#endif // MSL_MATRIX_VIEW_HPP