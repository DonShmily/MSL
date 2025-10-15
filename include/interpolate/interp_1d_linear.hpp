/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \include\interpolate\interp_1d_linear.hpp
** -----
** File Created: Wednesday, 15th October 2025 08:57:54
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Wednesday, 15th October 2025 14:08:29
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_INTERP_1D_LINEAR
#define MSL_INTERP_1D_LINEAR

#include "interp_1d.hpp"

namespace msl::interp
{
// Linear Interpolation
class Linear : public Interpolator1D
{
public:
    Linear() = default;

    Linear(std::span<const double> x, std::span<const double> y)
    {
        x_.assign(x.begin(), x.end());
        y_.assign(y.begin(), y.end());
        validate_input();
    }

    Linear(const std::vector<double> &x, const std::vector<double> &y)
        : Linear(std::span<const double>(x), std::span<const double>(y))
    {}

    void set_data(std::span<const double> x, std::span<const double> y)
    {
        x_.assign(x.begin(), x.end());
        y_.assign(y.begin(), y.end());
        validate_input();
    }

    double operator()(double x) const override
    {
        size_t i = find_interval(x);

        // Linear interpolation: y = y0 + (y1-y0)/(x1-x0) * (x-x0)
        double t = (x - x_[i]) / (x_[i + 1] - x_[i]);
        return y_[i] + t * (y_[i + 1] - y_[i]);
    }

    std::vector<double> operator()(std::span<const double> x_new) const override
    {
        std::vector<double> result(x_new.size());
        for (size_t i = 0; i < x_new.size(); ++i)
        {
            result[i] = (*this)(x_new[i]);
        }
        return result;
    }
};

inline std::vector<double> interp1_linear(std::span<const double> x,
                                          std::span<const double> y,
                                          std::span<const double> x_new)
{
    return Linear(x, y)(x_new);
}

} // namespace msl::interp

#endif // MSL_INTERP_1D_LINEAR