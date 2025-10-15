/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \include\matrix\complex_matrix_owmed.hpp
** -----
** File Created: Saturday, 11th October 2025 21:20:57
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Saturday, 11th October 2025 22:53:46
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_COMPLEX_MATRIX_OWNED_HPP
#define MSL_COMPLEX_MATRIX_OWNED_HPP

#include <algorithm>
#include <functional>
#include <initializer_list>
#include <vector>

#include "complex_matrix_base.hpp"

namespace msl
{
class complex_matrix_view; // Forward declaration

class complex_matrix_owned : public complex_matrix_base
{
private:
    using Base = complex_matrix_base;
    std::vector<std::complex<double>> storage_; // Actual data storage

    // Update base span when storage changes
    void update_span()
    {
        this->data_ = std::span<std::complex<double>>(storage_);
    }

public:
    // --- Constructors ---
    // Default constructor
    complex_matrix_owned() : Base() {}

    // Size constructor (zero-initialized)
    complex_matrix_owned(size_t rows, size_t cols)
        : Base(), storage_(rows * cols)
    {
        this->rows_ = rows;
        this->cols_ = cols;
        update_span();
    }

    // Size + value constructor
    complex_matrix_owned(size_t rows,
                         size_t cols,
                         const std::complex<double> &value)
        : Base(), storage_(rows * cols, value)
    {
        this->rows_ = rows;
        this->cols_ = cols;
        update_span();
    }

    // Construct from span (deep copy)
    complex_matrix_owned(size_t rows,
                         size_t cols,
                         std::span<const std::complex<double>> data)
        : Base(), storage_(data.begin(), data.end())
    {
        assert(data.size() == rows * cols);
        this->rows_ = rows;
        this->cols_ = cols;
        update_span();
    }

    // Construct from initializer list (row-major input)
    complex_matrix_owned(size_t rows,
                         size_t cols,
                         std::initializer_list<std::complex<double>> values)
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
    // --- Copy Constructors ---
    complex_matrix_owned(const complex_matrix_owned &other);
    explicit complex_matrix_owned(const complex_matrix_view &view);

    // --- Copy operations (deep copy) ---
    complex_matrix_owned &operator=(const complex_matrix_owned &other);
    complex_matrix_owned &operator=(const complex_matrix_view &view);

    // --- Move operations (no-throw guarantee) ---
    complex_matrix_owned(complex_matrix_owned &&other) noexcept;
    complex_matrix_owned &operator=(complex_matrix_owned &&other) noexcept;

    // --- Destructor ---
    ~complex_matrix_owned() = default;

    // --- Size management ---
    void resize(size_t rows,
                size_t cols,
                const std::complex<double> &value = std::complex<double>())
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
    void fill(const std::complex<double> &value)
    {
        std::fill(storage_.begin(), storage_.end(), value);
    }

    void fill_zeros()
    {
        std::fill(storage_.begin(), storage_.end(), std::complex<double>{});
    }

    // --- Row/Column extraction (returns new matrix) ---
    [[nodiscard]] complex_matrix_owned get_row(size_t i) const;
    [[nodiscard]] complex_matrix_owned get_column(size_t j) const;

    // --- Submatrix extraction ---
    [[nodiscard]] complex_matrix_owned submatrix(size_t row_start,
                                                 size_t row_end,
                                                 size_t col_start,
                                                 size_t col_end) const;

    // --- In-place operations ---
    complex_matrix_owned &operator+=(const complex_matrix_base &other);
    complex_matrix_owned &operator-=(const complex_matrix_base &other);
    complex_matrix_owned &operator*=(const std::complex<double> &scalar);
    complex_matrix_owned &operator/=(const std::complex<double> &scalar);

    // --- Apply function to each element ---
    void apply(std::function<std::complex<double>(std::complex<double>)> func);

    // --- Swap ---
    void swap(complex_matrix_owned &other) noexcept;

    // --- Factory methods ---
    // Create diagonal matrix from vector
    [[nodiscard]] static complex_matrix_owned
    diagonal(std::span<const std::complex<double>> diag);
};

typedef complex_matrix_owned matrixc;

// --- Free functions ---
void swap(complex_matrix_owned &a, complex_matrix_owned &b) noexcept;

// Binary operators (return new matrix)
[[nodiscard]] complex_matrix_owned operator+(const complex_matrix_owned &a,
                                             const complex_matrix_owned &b);
[[nodiscard]] complex_matrix_owned operator-(const complex_matrix_owned &a,
                                             const complex_matrix_owned &b);
[[nodiscard]] complex_matrix_owned
operator*(const complex_matrix_owned &mat, const std::complex<double> &scalar);
[[nodiscard]] complex_matrix_owned operator*(const std::complex<double> &scalar,
                                             const complex_matrix_owned &mat);
[[nodiscard]] complex_matrix_owned
operator/(const complex_matrix_owned &mat, const std::complex<double> &scalar);

} // namespace msl

#endif // MSL_complex_MATRIX_OWNED_HPP