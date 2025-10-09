/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \include\matrix\matrix.hpp
** -----
** File Created: Thursday, 9th October 2025 21:35:42
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Thursday, 9th October 2025 21:35:49
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

// Description: Concrete matrix class with owning semantics

#ifndef MSL_MATRIX_HPP
#define MSL_MATRIX_HPP

#include "matrix_base.hpp"
#include <cassert>
#include <span>
#include <vector>

namespace msl
{

template <typename T> class matrix : public matrix_base<T>
{
  private:
    std::vector<T> data_{}; // column-major storage

  public:
    using matrix_base<T>::rows_;
    using matrix_base<T>::cols_;

    // --- Constructors ---
    matrix() noexcept = default;

    matrix(size_t rows, size_t cols) : matrix_base<T>(rows, cols), data_(rows * cols)
    {
    }

    matrix(size_t rows, size_t cols, const T &value) : matrix_base<T>(rows, cols), data_(rows * cols, value)
    {
    }

    matrix(size_t rows, size_t cols, std::span<const T> data)
        : matrix_base<T>(rows, cols), data_(data.begin(), data.end())
    {
        assert(data_.size() == rows * cols);
    }

    // Copy / move default
    matrix(const matrix &) = default;
    matrix(matrix &&) noexcept = default;
    matrix &operator=(const matrix &) = default;
    matrix &operator=(matrix &&) noexcept = default;
    matrix &operator=(std::span<const T> data) noexcept
    {
        assert(data.size() == rows_ * cols_);
        std::copy(data.begin(), data.end(), data_.begin());
        return *this;
    }
    ~matrix() override = default;

    // --- Element access ---
    T &operator()(size_t i, size_t j) noexcept override
    {
        assert(i < rows_ && j < cols_);
        return data_[j * rows_ + i]; // column-major
    }

    const T &operator()(size_t i, size_t j) const noexcept override
    {
        assert(i < rows_ && j < cols_);
        return data_[j * rows_ + i];
    }

    // --- Column view ---
    std::span<T> get_column(size_t j) noexcept override
    {
        assert(j < cols_);
        return std::span<T>(data_.data() + j * rows_, rows_);
    }

    std::span<const T> get_column(size_t j) const noexcept override
    {
        assert(j < cols_);
        return std::span<const T>(data_.data() + j * rows_, rows_);
    }

    // --- Data access ---
    T *data() noexcept override
    {
        return data_.data();
    }
    const T *data() const noexcept override
    {
        return data_.data();
    }

    // --- Resize / fill helpers ---
    void resize(size_t rows, size_t cols)
    {
        rows_ = rows;
        cols_ = cols;
        data_.resize(rows * cols);
    }

    void fill(const T &value)
    {
        std::fill(data_.begin(), data_.end(), value);
    }

    // --- Consistency check ---
    bool is_consistent() const noexcept
    {
        return data_.size() == rows_ * cols_;
    }
};

using matrixd = matrix<double>;
using matrixc = matrix<std::complex<double>>;

} // namespace msl

#endif // MSL_MATRIX_HPP
