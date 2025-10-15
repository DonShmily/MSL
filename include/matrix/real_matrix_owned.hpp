/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \include\matrix\real_matrix_owned.hpp
** -----
** File Created: Saturday, 11th October 2025 21:20:14
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Saturday, 11th October 2025 22:42:59
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_REAL_MATRIX_OWNED_HPP
#define MSL_REAL_MATRIX_OWNED_HPP

#include <algorithm>
#include <initializer_list>
#include <vector>

#include "real_matrix_base.hpp"

namespace msl
{
class real_matrix_view; // Forward declaration

class real_matrix_owned : public real_matrix_base
{
private:
    using Base = real_matrix_base;
    std::vector<double> storage_; // Actual data storage

    // Update base span when storage changes
    void update_span() { this->data_ = std::span<double>(storage_); }

public:
    // --- Constructors ---
    // Default constructor
    real_matrix_owned() : Base() {}

    // Size constructor (zero-initialized)
    real_matrix_owned(size_t rows, size_t cols) : Base(), storage_(rows * cols)
    {
        this->rows_ = rows;
        this->cols_ = cols;
        update_span();
    }

    // Size + value constructor
    real_matrix_owned(size_t rows, size_t cols, const double &value)
        : Base(), storage_(rows * cols, value)
    {
        this->rows_ = rows;
        this->cols_ = cols;
        update_span();
    }

    // Construct from span (deep copy)
    real_matrix_owned(size_t rows, size_t cols, std::span<const double> data)
        : Base(), storage_(data.begin(), data.end())
    {
        assert(data.size() == rows * cols);
        this->rows_ = rows;
        this->cols_ = cols;
        update_span();
    }

    // Construct from initializer list (row-major input)
    real_matrix_owned(size_t rows,
                      size_t cols,
                      std::initializer_list<double> values)
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
    real_matrix_owned(const real_matrix_owned &other);
    explicit real_matrix_owned(const real_matrix_view &view);

    // --- Copy operations (deep copy) ---
    real_matrix_owned &operator=(const real_matrix_owned &other);
    real_matrix_owned &operator=(const real_matrix_view &view);

    // --- Move operations (no-throw guarantee) ---
    real_matrix_owned(real_matrix_owned &&other) noexcept;
    real_matrix_owned &operator=(real_matrix_owned &&other) noexcept;

    // --- Destructor ---
    ~real_matrix_owned() = default;

    // --- Size management ---
    void resize(size_t rows, size_t cols, const double &value = double())
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
    void fill(const double &value)
    {
        std::fill(storage_.begin(), storage_.end(), value);
    }

    void fill_zeros() { std::fill(storage_.begin(), storage_.end(), double{}); }

    // --- Row/Column extraction (returns new matrix) ---
    [[nodiscard]] real_matrix_owned get_row(size_t i) const;
    [[nodiscard]] real_matrix_owned get_column(size_t j) const;

    // --- Submatrix extraction ---
    [[nodiscard]] real_matrix_owned submatrix(size_t row_start,
                                              size_t row_end,
                                              size_t col_start,
                                              size_t col_end) const;

    // --- In-place operations ---
    real_matrix_owned &operator+=(const real_matrix_base &other);
    real_matrix_owned &operator-=(const real_matrix_base &other);
    real_matrix_owned &operator*=(const double &scalar);
    real_matrix_owned &operator/=(const double &scalar);

    // --- Swap ---
    void swap(real_matrix_owned &other) noexcept;

    // --- Factory methods ---
    [[nodiscard]] static real_matrix_owned zeros(size_t rows, size_t cols);
    [[nodiscard]] static real_matrix_owned ones(size_t rows, size_t cols);
    [[nodiscard]] static real_matrix_owned identity(size_t n);
    // Create diagonal matrix from vector
    [[nodiscard]] static real_matrix_owned
    diagonal(std::span<const double> diag);
};

typedef real_matrix_owned matrixd;

// --- Free functions ---
void swap(real_matrix_owned &a, real_matrix_owned &b) noexcept;

// Binary operators (return new matrix)
[[nodiscard]] real_matrix_owned operator+(const real_matrix_owned &a,
                                          const real_matrix_owned &b);
[[nodiscard]] real_matrix_owned operator-(const real_matrix_owned &a,
                                          const real_matrix_owned &b);
[[nodiscard]] real_matrix_owned operator*(const real_matrix_owned &mat,
                                          const double &scalar);
[[nodiscard]] real_matrix_owned operator*(const double &scalar,
                                          const real_matrix_owned &mat);
[[nodiscard]] real_matrix_owned operator/(const real_matrix_owned &mat,
                                          const double &scalar);

} // namespace msl

#endif // MSL_REAL_MATRIX_OWNED_HPP