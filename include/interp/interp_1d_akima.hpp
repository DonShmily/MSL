/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: interp_1d_akima.hpp
** -----
** File Created: Wednesday, 15th October 2025 14:12:29
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:03:33
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_INTERP_1D_AKIMA
#define MSL_INTERP_1D_AKIMA

#include "interp_1d_base.hpp"

namespace msl::interp
{
// Akima Spline (smoother than cubic, no overshoot)
class AkimaSpline : public InterpolatorBase
{
private:
    std::vector<double> b_, c_, d_;

    void compute_coefficients()
    {
        size_t n = x_.size();

        // Compute slopes
        std::vector<double> m(n + 3);

        // Interior slopes
        for (size_t i = 2; i < n + 1; ++i)
        {
            m[i] = (y_[i - 1] - y_[i - 2]) / (x_[i - 1] - x_[i - 2]);
        }

        // Extrapolate for boundary
        m[0] = 3.0 * m[2] - 2.0 * m[3];
        m[1] = 2.0 * m[2] - m[3];
        m[n + 1] = 2.0 * m[n] - m[n - 1];
        m[n + 2] = 3.0 * m[n] - 2.0 * m[n - 1];

        // Compute Akima weights
        std::vector<double> t(n);
        for (size_t i = 0; i < n; ++i)
        {
            double w1 = std::abs(m[i + 3] - m[i + 2]);
            double w2 = std::abs(m[i + 1] - m[i]);

            if (w1 + w2 < 1e-10)
            {
                t[i] = 0.5 * (m[i + 1] + m[i + 2]);
            }
            else
            {
                t[i] = (w1 * m[i + 1] + w2 * m[i + 2]) / (w1 + w2);
            }
        }

        // Compute polynomial coefficients
        b_.resize(n - 1);
        c_.resize(n - 1);
        d_.resize(n - 1);

        for (size_t i = 0; i < n - 1; ++i)
        {
            double h = x_[i + 1] - x_[i];
            b_[i] = t[i];
            c_[i] = (3.0 * m[i + 2] - 2.0 * t[i] - t[i + 1]) / h;
            d_[i] = (t[i] + t[i + 1] - 2.0 * m[i + 2]) / (h * h);
        }
    }

public:
    AkimaSpline() = default;

    AkimaSpline(std::span<const double> x, std::span<const double> y)
    {
        x_.assign(x.begin(), x.end());
        y_.assign(y.begin(), y.end());
        validate_input();
        if (x_.size() < 4)
        {
            throw std::invalid_argument("Akima spline needs at least 4 points");
        }
        compute_coefficients();
    }

    void set_data(std::span<const double> x, std::span<const double> y) override
    {
        x_.assign(x.begin(), x.end());
        y_.assign(y.begin(), y.end());
        validate_input();
        if (x_.size() < 4)
        {
            throw std::invalid_argument("Akima spline needs at least 4 points");
        }
        compute_coefficients();
    }

    double interpolate(double x) const override
    {
        size_t i = find_interval(x);
        double dx = x - x_[i];
        return y_[i] + b_[i] * dx + c_[i] * dx * dx + d_[i] * dx * dx * dx;
    }
};

inline std::vector<double> interp1_akima(std::span<const double> x,
                                         std::span<const double> y,
                                         std::span<const double> x_new)
{
    return AkimaSpline(x, y)(x_new);
}

} // namespace msl::interp
#endif // MSL_INTERP_1D_AKIMA