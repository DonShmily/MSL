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
** Last Modified: Thursday, 9th October 2025 22:30:18
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

// Description: Non-owning matrix view class

#ifndef MSL_MATRIX_VIEW_HPP
#define MSL_MATRIX_VIEW_HPP

#include "matrix_base.hpp"
#include <cassert>
#include <span>

namespace msl
{

template <typename T> class matrix_view : public matrix_base<T>
{
  private:
    T *ptr_{nullptr};

  public:
    using matrix_base<T>::rows_;
    using matrix_base<T>::cols_;

    // --- Constructors ---
    matrix_view() noexcept = default;

    matrix_view(T *ptr, size_t rows, size_t cols) noexcept : matrix_base<T>(rows, cols), ptr_(ptr)
    {
    }

    // Construct from span (for convenience)
    matrix_view(std::span<T> data, size_t rows, size_t cols) noexcept : matrix_base<T>(rows, cols), ptr_(data.data())
    {
        assert(data.size() == rows * cols);
    }

    // --- Element access ---
    T &operator()(size_t i, size_t j) noexcept override
    {
        assert(i < rows_ && j < cols_);
        return ptr_[j * rows_ + i];
    }

    const T &operator()(size_t i, size_t j) const noexcept override
    {
        assert(i < rows_ && j < cols_);
        return ptr_[j * rows_ + i];
    }

    // --- Column view ---
    std::span<T> get_column(size_t j) noexcept override
    {
        assert(j < cols_);
        return std::span<T>(ptr_ + j * rows_, rows_);
    }

    std::span<const T> get_column(size_t j) const noexcept override
    {
        assert(j < cols_);
        return std::span<const T>(ptr_ + j * rows_, rows_);
    }

    // --- Data access ---
    T *data() noexcept override
    {
        return ptr_;
    }
    const T *data() const noexcept override
    {
        return ptr_;
    }

    // --- Validity check ---
    bool valid() const noexcept
    {
        return ptr_ != nullptr && rows_ > 0 && cols_ > 0;
    }
};

} // namespace msl

#endif // MSL_MATRIX_VIEW_HPP
