/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025
**  File: interpolator_base.hpp
**  Description: Base class for interpolators with extrapolation support
*/

#ifndef MSL_INTERPOLATOR_BASE_HPP
#define MSL_INTERPOLATOR_BASE_HPP

#include <algorithm>
#include <cmath>
#include <span>
#include <stdexcept>
#include <vector>


namespace msl::interp
{

/**
 * @brief Extrapolation methods
 */
enum class ExtrapolationMode
{
    None,      // Throw exception (default)
    Constant,  // Use boundary values
    Linear,    // Linear extrapolation using boundary slopes
    Nearest,   // Use nearest boundary point
    Periodic,  // Periodic extension
    Polynomial // Use boundary polynomial (MATLAB-style for spline/Akima)
};

/**
 * @brief Base class for 1D interpolators
 */
class InterpolatorBase
{
protected:
    std::vector<double> x_;
    std::vector<double> y_;
    ExtrapolationMode extrap_mode_ = ExtrapolationMode::Polynomial;

    void validate_input() const
    {
        if (x_.size() != y_.size())
        {
            throw std::invalid_argument("x and y must have same size");
        }
        if (x_.size() < 2)
        {
            throw std::invalid_argument("Need at least 2 points");
        }

        // Check sorted
        for (size_t i = 1; i < x_.size(); ++i)
        {
            if (x_[i] <= x_[i - 1])
            {
                throw std::invalid_argument(
                    "x values must be strictly increasing");
            }
        }
    }

    /**
     * @brief Find interval for interpolation point
     * Throws if out of range and extrapolation is disabled
     */
    size_t find_interval(double x) const
    {
        if (x < x_.front())
        {
            if (extrap_mode_ == ExtrapolationMode::None)
            {
                throw std::out_of_range("x out of interpolation range");
            }
            // Handle extrapolation
            if (extrap_mode_ == ExtrapolationMode::Polynomial)
                return 0;
        }
        if (x > x_.back())
        {
            if (extrap_mode_ == ExtrapolationMode::None)
            {
                throw std::out_of_range("x out of interpolation range");
            }
            // Handle extrapolation
            if (extrap_mode_ == ExtrapolationMode::Polynomial)
                return x_.size() - 2;
        }

        // Binary search
        auto it = std::lower_bound(x_.begin(), x_.end(), x);
        size_t idx = it - x_.begin();

        if (idx == x_.size())
            --idx;
        if (idx > 0 && x < x_[idx])
            --idx;

        return idx;
    }

    /**
     * @brief Pure virtual interpolation function
     * Must be implemented by derived classes
     */
    virtual double interpolate(double x) const = 0;

public:
    InterpolatorBase() = default;
    virtual ~InterpolatorBase() = default;

    /**
     * @brief Set extrapolation mode
     */
    void set_extrapolation(ExtrapolationMode mode) { extrap_mode_ = mode; }

    /**
     * @brief Get current extrapolation mode
     */
    [[nodiscard]] ExtrapolationMode extrapolation_mode() const
    {
        return extrap_mode_;
    }

    /**
     * @brief Single point evaluation
     */
    double operator()(double x) const
    {
        // Check for exact match
        auto idx = std::find(x_.begin(), x_.end(), x);
        if (idx != x_.end())
        {
            return y_[idx - x_.begin()];
        }

        // Handle extrapolation
        if (x < x_.front())
        {
            switch (extrap_mode_)
            {
                case ExtrapolationMode::None:
                    throw std::out_of_range(
                        "x < x_min, extrapolation disabled");

                case ExtrapolationMode::Constant:
                    return y_.front();

                case ExtrapolationMode::Linear: {
                    // Linear extrapolation using first segment slope
                    double slope = (y_[1] - y_[0]) / (x_[1] - x_[0]);
                    return y_[0] + slope * (x - x_[0]);
                }

                case ExtrapolationMode::Nearest:
                    return y_.front();

                case ExtrapolationMode::Periodic: {
                    double period = x_.back() - x_.front();
                    double x_wrapped =
                        x + period * std::ceil((x_.front() - x) / period);
                    return interpolate(x_wrapped);
                }

                case ExtrapolationMode::Polynomial: {
                    break;
                    ; // Let derived class handle polynomial extrapolation
                }
            }
        }

        if (x > x_.back())
        {
            switch (extrap_mode_)
            {
                case ExtrapolationMode::None:
                    throw std::out_of_range(
                        "x > x_max, extrapolation disabled");

                case ExtrapolationMode::Constant:
                    return y_.back();

                case ExtrapolationMode::Linear: {
                    // Linear extrapolation using last segment slope
                    size_t n = x_.size();
                    double slope =
                        (y_[n - 1] - y_[n - 2]) / (x_[n - 1] - x_[n - 2]);
                    return y_[n - 1] + slope * (x - x_[n - 1]);
                }

                case ExtrapolationMode::Nearest:
                    return y_.back();

                case ExtrapolationMode::Periodic: {
                    double period = x_.back() - x_.front();
                    double x_wrapped =
                        x - period * std::floor((x - x_.front()) / period);
                    return interpolate(x_wrapped);
                }

                case ExtrapolationMode::Polynomial: {
                    break;
                    ; // Let derived class handle polynomial extrapolation
                }
            }
        }

        // Within range - use derived class interpolation
        return interpolate(x);
    }

    /**
     * @brief Multiple point evaluation
     */
    std::vector<double> operator()(std::span<const double> x_new) const
    {
        std::vector<double> result(x_new.size());
        for (size_t i = 0; i < x_new.size(); ++i)
        {
            result[i] = (*this)(x_new[i]);
        }
        return result;
    }

    std::vector<double> operator()(const std::vector<double> &x_new) const
    {
        return (*this)(std::span<const double>(x_new));
    }

    /**
     * @brief Set data points (pure virtual)
     */
    virtual void set_data(std::span<const double> x,
                          std::span<const double> y) = 0;

    /**
     * @brief Get data points
     */
    [[nodiscard]] const std::vector<double> &x_data() const { return x_; }
    [[nodiscard]] const std::vector<double> &y_data() const { return y_; }
    [[nodiscard]] size_t size() const { return x_.size(); }

    /**
     * @brief Check if point is within interpolation range
     */
    [[nodiscard]] bool in_range(double x) const
    {
        return x >= x_.front() && x <= x_.back();
    }

    /**
     * @brief Get interpolation range
     */
    [[nodiscard]] std::pair<double, double> range() const
    {
        return {x_.front(), x_.back()};
    }
};

} // namespace msl::interp

#endif // MSL_INTERPOLATOR_BASE_HPP