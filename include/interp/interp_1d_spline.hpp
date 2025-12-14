/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: interp_1d_spline.hpp
** -----
** File Created: Wednesday, 15th October 2025 14:11:26
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:03:56
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_INTERP_1D_SPLINE
#define MSL_INTERP_1D_SPLINE

#include "interp_1d_base.hpp"

namespace msl::interp
{
// Cubic Spline Interpolation
class CubicSpline : public InterpolatorBase
{
private:
    std::vector<double> b_; // First derivative coefficients
    std::vector<double> c_; // Second derivative coefficients
    std::vector<double> d_; // Third derivative coefficients

    void compute_coefficients()
    {
        size_t n = x_.size();

        // Compute intervals
        std::vector<double> h(n - 1);
        for (size_t i = 0; i < n - 1; ++i)
        {
            h[i] = x_[i + 1] - x_[i];
        }

        // Build tridiagonal system for natural spline
        // Natural boundary: c[0] = c[n-1] = 0
        std::vector<double> alpha(n);
        for (size_t i = 1; i < n - 1; ++i)
        {
            alpha[i] =
                3.0
                * ((y_[i + 1] - y_[i]) / h[i] - (y_[i] - y_[i - 1]) / h[i - 1]);
        }

        // Solve tridiagonal system using Thomas algorithm
        std::vector<double> l(n), mu(n), z(n);
        l[0] = 1.0;
        mu[0] = 0.0;
        z[0] = 0.0;

        for (size_t i = 1; i < n - 1; ++i)
        {
            l[i] = 2.0 * (x_[i + 1] - x_[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }

        l[n - 1] = 1.0;
        z[n - 1] = 0.0;

        // Back substitution
        c_.resize(n);
        c_[n - 1] = 0.0;

        for (int j = n - 2; j >= 0; --j)
        {
            c_[j] = z[j] - mu[j] * c_[j + 1];
        }

        // Compute b and d coefficients
        b_.resize(n - 1);
        d_.resize(n - 1);

        for (size_t i = 0; i < n - 1; ++i)
        {
            b_[i] = (y_[i + 1] - y_[i]) / h[i]
                    - h[i] * (c_[i + 1] + 2.0 * c_[i]) / 3.0;
            d_[i] = (c_[i + 1] - c_[i]) / (3.0 * h[i]);
        }
    }

public:
    CubicSpline() = default;

    CubicSpline(std::span<const double> x, std::span<const double> y)
    {
        x_.assign(x.begin(), x.end());
        y_.assign(y.begin(), y.end());
        validate_input();
        compute_coefficients();
    }

    CubicSpline(const std::vector<double> &x, const std::vector<double> &y)
        : CubicSpline(std::span<const double>(x), std::span<const double>(y))
    {}

    void set_data(std::span<const double> x, std::span<const double> y) override
    {
        x_.assign(x.begin(), x.end());
        y_.assign(y.begin(), y.end());
        validate_input();
        compute_coefficients();
    }

    double interpolate(double x) const override
    {
        size_t i = find_interval(x);

        // Cubic spline: S_i(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 +
        // d_i*(x-x_i)^3
        double dx = x - x_[i];
        return y_[i] + b_[i] * dx + c_[i] * dx * dx + d_[i] * dx * dx * dx;
    }

    // Get derivative at point
    double derivative(double x) const
    {
        size_t i = find_interval(x);
        double dx = x - x_[i];
        return b_[i] + 2.0 * c_[i] * dx + 3.0 * d_[i] * dx * dx;
    }
};

inline std::vector<double> interp1_cubic(std::span<const double> x,
                                         std::span<const double> y,
                                         std::span<const double> x_new)
{
    return CubicSpline(x, y)(x_new);
}

} // namespace msl::interp

#endif // MSL_INTERP_1D_SPLINE