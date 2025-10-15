/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \include\matrix\matrix_base.hpp
** -----
** File Created: Thursday, 9th October 2025 09:55:48
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Friday, 10th October 2025 17:41:24
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/


#ifndef MSL_MATRIX_BASE_HPP
#define MSL_MATRIX_BASE_HPP

#include <cassert>
#include <span>
#include <utility>


// Optional bounds checking (enable with -DMSL_BOUNDS_CHECK)
#ifdef MSL_BOUNDS_CHECK
#define MSL_CHECK_BOUNDS(condition, message)                                   \
    if (!(condition)) throw std::out_of_range(message)
#else
#define MSL_CHECK_BOUNDS(condition, message) assert(condition)
#endif

namespace msl
{

/**
 * @brief CRTP base class for matrix types
 *
 * Uses Curiously Recurring Template Pattern to eliminate virtual function
 * overhead while maintaining a unified interface.
 *
 * @tparam T Element type (double, std::complex<double>, etc.)
 * @tparam Derived The derived class (matrix_owned<T> or matrix_view<T>)
 */
template <typename T, typename Derived>
class matrix_base
{
protected:
    std::span<T> data_; // Unified data access through span
    size_t rows_{0};
    size_t cols_{0};

    // Protected constructor - only derived classes can construct
    matrix_base() = default;

    matrix_base(std::span<T> data, size_t rows, size_t cols)
        : data_(data), rows_(rows), cols_(cols)
    {
        assert(data_.size() == rows * cols);
    }

public:
    // --- Type traits ---
    using value_type = T;
    using reference = T &;
    using const_reference = const T &;
    using pointer = T *;
    using const_pointer = const T *;
    using size_type = size_t;

    // --- Size queries ---
    [[nodiscard]] size_t rows() const noexcept { return rows_; }
    [[nodiscard]] size_t cols() const noexcept { return cols_; }
    [[nodiscard]] size_t size() const noexcept { return rows_ * cols_; }
    [[nodiscard]] std::pair<size_t, size_t> shape() const noexcept
    {
        return {rows_, cols_};
    }
    [[nodiscard]] bool empty() const noexcept { return size() == 0; }

    // --- Element access (CRTP: no virtual function overhead) ---
    [[nodiscard]] T &operator()(size_t i, size_t j) noexcept
    {
        MSL_CHECK_BOUNDS(i < rows_ && j < cols_, "Matrix index out of bounds");
        return data_[j * rows_ + i]; // Column-major order
    }

    [[nodiscard]] const T &operator()(size_t i, size_t j) const noexcept
    {
        MSL_CHECK_BOUNDS(i < rows_ && j < cols_, "Matrix index out of bounds");
        return data_[j * rows_ + i];
    }

    // Linear indexing (useful for iteration)
    [[nodiscard]] T &operator[](size_t idx) noexcept
    {
        MSL_CHECK_BOUNDS(idx < size(), "Linear index out of bounds");
        return data_[idx];
    }

    [[nodiscard]] const T &operator[](size_t idx) const noexcept
    {
        MSL_CHECK_BOUNDS(idx < size(), "Linear index out of bounds");
        return data_[idx];
    }

    // --- Column access (efficient in column-major layout) ---
    [[nodiscard]] std::span<T> column(size_t j) noexcept
    {
        MSL_CHECK_BOUNDS(j < cols_, "Column index out of bounds");
        return {data_.data() + j * rows_, rows_};
    }

    [[nodiscard]] std::span<const T> column(size_t j) const noexcept
    {
        MSL_CHECK_BOUNDS(j < cols_, "Column index out of bounds");
        return {data_.data() + j * rows_, rows_};
    }

    // --- Data access ---
    [[nodiscard]] T *data() noexcept { return data_.data(); }
    [[nodiscard]] const T *data() const noexcept { return data_.data(); }
    [[nodiscard]] std::span<T> data_span() noexcept { return data_; }
    [[nodiscard]] std::span<const T> data_span() const noexcept
    {
        return data_;
    }

    // --- Iterators ---
    [[nodiscard]] auto begin() noexcept { return data_.begin(); }
    [[nodiscard]] auto end() noexcept { return data_.end(); }
    [[nodiscard]] auto begin() const noexcept { return data_.begin(); }
    [[nodiscard]] auto end() const noexcept { return data_.end(); }
    [[nodiscard]] auto cbegin() const noexcept { return data_.cbegin(); }
    [[nodiscard]] auto cend() const noexcept { return data_.cend(); }


    // --- Derived class access (CRTP) ---
    [[nodiscard]] Derived &derived() noexcept
    {
        return static_cast<Derived &>(*this);
    }

    [[nodiscard]] const Derived &derived() const noexcept
    {
        return static_cast<const Derived &>(*this);
    }

protected:
    // Helper for checking consistency
    [[nodiscard]] bool is_consistent() const noexcept
    {
        return data_.size() == rows_ * cols_;
    }
};

} // namespace msl

#endif // MSL_MATRIX_BASE_HPP