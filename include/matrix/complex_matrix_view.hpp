/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \include\matrix\complex_matrix_view.hpp
** -----
** File Created: Saturday, 11th October 2025 21:20:53
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Saturday, 11th October 2025 22:49:24
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/


#ifndef MSL_complex_MATRIX_VIEW_HPP
#define MSL_complex_MATRIX_VIEW_HPP

#include <cassert>
#include "complex_matrix_base.hpp"
#include "complex_matrix_owned.hpp"

namespace msl
{
class const_complex_matrix_view; // Forward declaration

class complex_matrix_view : public complex_matrix_base
{
private:
    using Base = complex_matrix_base;

public:
    // --- Constructors ---
    // Default constructor (empty view)
    complex_matrix_view() : Base() {}

    // Construct from raw pointer
    complex_matrix_view(std::complex<double> *data, size_t rows, size_t cols)
        : Base(rows, cols, std::span<std::complex<double>>(data, rows * cols))
    {}

    // Construct from complex_matrix_owned (most common use case)
    complex_matrix_view(complex_matrix_owned &owner)
        : Base(owner.rows(), owner.cols(), owner.data_span())
    {}

    // Construct from span
    complex_matrix_view(std::span<std::complex<double>> data,
                        size_t rows,
                        size_t cols)
        : Base(rows, cols, data)
    {
        assert(data.size() == rows * cols);
    }

    // --- Copy semantics: Shallow copy (copies the view) ---
    complex_matrix_view(const complex_matrix_view &other) = default;
    complex_matrix_view &operator=(const complex_matrix_view &other) = default;

    // --- Move semantics ---
    complex_matrix_view(complex_matrix_view &&other) noexcept = default;
    complex_matrix_view &
    operator=(complex_matrix_view &&other) noexcept = default;

    // --- Cannot construct view from const owned (use const_matrix_view) ---
    complex_matrix_view(const complex_matrix_owned &) = delete;

    // --- Destructor ---
    ~complex_matrix_view() = default;

    // --- Subview creation ---
    [[nodiscard]] complex_matrix_view
    subview(size_t row_start, size_t row_end, size_t col_start, size_t col_end);

    // --- Column subview (efficient) ---
    [[nodiscard]] complex_matrix_view column_view(size_t j);
    [[nodiscard]] const_complex_matrix_view column_view(size_t j) const;

    // --- Fill operations ---
    void fill(const std::complex<double> &value)
    {
        std::fill(this->begin(), this->end(), value);
    }

    void fill_zeros()
    {
        std::fill(this->begin(), this->end(), std::complex<double>{0.0, 0.0});
    }

    // --- In-place operations ---
    complex_matrix_view &operator+=(const complex_matrix_base &other);
    complex_matrix_view &operator-=(const complex_matrix_base &other);
    complex_matrix_view &operator*=(const std::complex<double> &scalar);
    complex_matrix_view &operator/=(const std::complex<double> &scalar);

    // --- Copy data from another matrix ---
    void copy_from(const complex_matrix_base &src)
    {
        assert(this->rows_ == src.rows() && this->cols_ == src.cols());
        std::copy(src.begin(), src.end(), this->begin());
    }

    // --- Convert to owned matrix (deep copy) ---
    [[nodiscard]] complex_matrix_owned to_owned() const
    {
        return complex_matrix_owned(this->rows_, this->cols_, this->data_);
    }

    // --- Swap (swaps the views, not the data) ---
    void swap(complex_matrix_view &other) noexcept
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
class const_complex_matrix_view
{
private:
    std::span<const std::complex<double>> data_;
    size_t rows_{0};
    size_t cols_{0};

public:
    // --- Constructors ---
    const_complex_matrix_view() = default;

    const_complex_matrix_view(const std::complex<double> *data,
                              size_t rows,
                              size_t cols)
        : data_(data, rows * cols), rows_(rows), cols_(cols)
    {}

    const_complex_matrix_view(const complex_matrix_owned &owner)
        : data_(owner.data_span()), rows_(owner.rows()), cols_(owner.cols())
    {}

    const_complex_matrix_view(const complex_matrix_view &view)
        : data_(view.data_span()), rows_(view.rows()), cols_(view.cols())
    {}

    const_complex_matrix_view(std::span<const std::complex<double>> data,
                              size_t rows,
                              size_t cols)
        : data_(data), rows_(rows), cols_(cols)
    {
        assert(data.size() == rows * cols);
    }

    // --- Copy and move ---
    const_complex_matrix_view(const const_complex_matrix_view &) = default;
    const_complex_matrix_view &
    operator=(const const_complex_matrix_view &) = default;
    const_complex_matrix_view(const_complex_matrix_view &&) noexcept = default;
    const_complex_matrix_view &
    operator=(const_complex_matrix_view &&) noexcept = default;

    // --- Size queries ---
    [[nodiscard]] size_t rows() const noexcept { return rows_; }
    [[nodiscard]] size_t cols() const noexcept { return cols_; }
    [[nodiscard]] size_t size() const noexcept { return rows_ * cols_; }
    [[nodiscard]] std::pair<size_t, size_t> shape() const noexcept
    {
        return {rows_, cols_};
    }

    // --- Element access (read-only) ---
    [[nodiscard]] const std::complex<double> &
    operator()(size_t i, size_t j) const noexcept
    {
        return data_[j * rows_ + i];
    }

    [[nodiscard]] const std::complex<double> &
    operator[](size_t idx) const noexcept
    {
        return data_[idx];
    }

    // --- Column access ---
    [[nodiscard]] std::span<const std::complex<double>>
    column(size_t j) const noexcept
    {
        return {data_.data() + j * rows_, rows_};
    }

    // --- Data access ---
    [[nodiscard]] const std::complex<double> *data() const noexcept
    {
        return data_.data();
    }
    [[nodiscard]] std::span<const std::complex<double>>
    data_span() const noexcept
    {
        return data_;
    }

    // --- Iterators ---
    [[nodiscard]] auto begin() const noexcept { return data_.begin(); }
    [[nodiscard]] auto end() const noexcept { return data_.end(); }
    [[nodiscard]] auto cbegin() const noexcept { return data_.begin(); }
    [[nodiscard]] auto cend() const noexcept { return data_.end(); }

    // --- Convert to owned matrix ---
    [[nodiscard]] complex_matrix_owned to_owned() const
    {
        return complex_matrix_owned(
            this->rows_,
            this->cols_,
            std::span<const std::complex<double>>(this->data_));
    }
};

typedef complex_matrix_view matrixc_view;
typedef const_complex_matrix_view const_matrixc_view;

// --- Helper functions ---
void swap(complex_matrix_view &a, complex_matrix_view &b) noexcept;
// Create view from owned matrix
[[nodiscard]] complex_matrix_view make_view(complex_matrix_owned &mat);
[[nodiscard]] const_complex_matrix_view
make_view(const complex_matrix_owned &mat);

} // namespace msl

#endif // MSL_MATRIX_VIEW_HPP