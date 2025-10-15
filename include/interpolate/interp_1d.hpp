/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025
**  File: interpolation.hpp
**  Description: 1D interpolation (linear, cubic spline)
*/

#ifndef MSL_INTERP_1D
#define MSL_INTERP_1D

#include <algorithm>
#include <span>
#include <stdexcept>
#include <vector>


namespace msl::interp
{
// Base class for 1D interpolation
class Interpolator1D
{
protected:
    std::vector<double> x_;
    std::vector<double> y_;

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

    // Binary search for interval
    size_t find_interval(double x) const
    {
        if (x < x_.front() || x > x_.back())
        {
            throw std::out_of_range("x out of interpolation range");
        }

        // Binary search
        auto it = std::lower_bound(x_.begin(), x_.end(), x);
        size_t idx = it - x_.begin();

        if (idx == x_.size())
            --idx;
        if (idx > 0 && (x < x_[idx]))
            --idx;

        return idx;
    }

public:
    Interpolator1D() = default;
    virtual ~Interpolator1D() = default;

    virtual double operator()(double x) const = 0;
    virtual std::vector<double>
    operator()(std::span<const double> x_new) const = 0;

    [[nodiscard]] const std::vector<double> &x_data() const { return x_; }
    [[nodiscard]] const std::vector<double> &y_data() const { return y_; }
    [[nodiscard]] size_t size() const { return x_.size(); }
};
} // namespace msl::interp

#endif // MSL_INTERP_1D