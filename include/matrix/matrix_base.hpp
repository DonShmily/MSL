// Description: Abstract base class for matrix types (owning and non-owning)

#ifndef MSL_MATRIX_BASE_HPP
#define MSL_MATRIX_BASE_HPP

#include <cstddef>
#include <span>
#include <utility>

namespace msl
{

// Base class defining a common interface for matrix-like objects
template <typename T> class matrix_base
{
  protected:
    size_t rows_{0}, cols_{0}; // dimensions

  public:
    // Default constructor (empty matrix)
    constexpr matrix_base() noexcept = default;

    // Constructor with specified shape
    constexpr matrix_base(size_t rows, size_t cols) noexcept : rows_(rows), cols_(cols)
    {
    }

    // Virtual destructor (required for base class)
    virtual ~matrix_base() = default;

    // --- Basic info ---
    [[nodiscard]] constexpr size_t rows() const noexcept
    {
        return rows_;
    }
    [[nodiscard]] constexpr size_t cols() const noexcept
    {
        return cols_;
    }

    [[nodiscard]] constexpr std::pair<size_t, size_t> size() const noexcept
    {
        return {rows_, cols_};
    }

    [[nodiscard]] constexpr bool empty() const noexcept
    {
        return rows_ == 0 || cols_ == 0;
    }

    // --- Element access interface (to be overridden) ---
    virtual T &operator()(size_t i, size_t j) noexcept = 0;
    virtual const T &operator()(size_t i, size_t j) const noexcept = 0;

    // --- Column view interface (to be overridden) ---
    virtual std::span<T> get_column(size_t j) noexcept = 0;
    virtual std::span<const T> get_column(size_t j) const noexcept = 0;

    // --- Raw data access (to be overridden) ---
    virtual T *data() noexcept = 0;
    virtual const T *data() const noexcept = 0;
};

} // namespace msl

#endif // MSL_MATRIX_BASE_HPP
