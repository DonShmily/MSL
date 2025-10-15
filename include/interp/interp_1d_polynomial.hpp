/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025
**  File: polynomial_interpolator.hpp
**  Description: Polynomial interpolation (Lagrange & Newton)
*/

#ifndef MSL_POLYNOMIAL_INTERPOLATOR_HPP
#define MSL_POLYNOMIAL_INTERPOLATOR_HPP

#include <span>
#include <vector>

#include "interp_1d_base.hpp"

namespace msl::interp
{

/**
 * @brief Polynomial interpolation using divided differences (Newton form)
 *
 * More stable than Lagrange form, allows incremental updates.
 * Warning: High degree polynomials can oscillate (Runge phenomenon)
 *
 * Recommended: Use for <= 10 points, or consider splines instead
 */
class PolynomialInterpolator : public InterpolatorBase
{
private:
    std::vector<double> coeffs_; // Divided difference coefficients

    void compute_divided_differences()
    {
        size_t n = x_.size();
        coeffs_.resize(n);

        // Initialize with y values
        std::vector<std::vector<double>> table(n);
        for (size_t i = 0; i < n; ++i)
        {
            table[i].resize(n - i);
            table[i][0] = y_[i];
        }

        // Compute divided differences
        for (size_t j = 1; j < n; ++j)
        {
            for (size_t i = 0; i < n - j; ++i)
            {
                table[i][j] = (table[i + 1][j - 1] - table[i][j - 1])
                              / (x_[i + j] - x_[i]);
            }
        }

        // Extract coefficients (first row)
        for (size_t i = 0; i < n; ++i)
        {
            coeffs_[i] = table[0][i];
        }
    }

public:
    PolynomialInterpolator() = default;

    PolynomialInterpolator(std::span<const double> x, std::span<const double> y)
    {
        x_.assign(x.begin(), x.end());
        y_.assign(y.begin(), y.end());
        validate_input();
        compute_divided_differences();
    }

    PolynomialInterpolator(const std::vector<double> &x,
                           const std::vector<double> &y)
        : PolynomialInterpolator(std::span<const double>(x),
                                 std::span<const double>(y))
    {}

    void set_data(std::span<const double> x, std::span<const double> y) override
    {
        x_.assign(x.begin(), x.end());
        y_.assign(y.begin(), y.end());
        validate_input();
        compute_divided_differences();
    }

    double interpolate(double x) const override
    {
        // Newton's form: P(x) = c0 + c1(x-x0) + c2(x-x0)(x-x1) + ...
        size_t n = x_.size();
        double result = coeffs_[n - 1];

        // Horner's method (reverse order for stability)
        for (int i = n - 2; i >= 0; --i)
        {
            result = result * (x - x_[i]) + coeffs_[i];
        }

        return result;
    }

    /**
     * @brief Get polynomial degree
     */
    [[nodiscard]] size_t degree() const { return x_.size() - 1; }

    /**
     * @brief Evaluate derivative at point
     *
     * Computes derivative using direct differentiation of Newton form
     */
    double derivative(double x) const
    {
        size_t n = x_.size();
        if (n < 2)
            return 0.0;

        double result = 0.0;

        for (size_t k = 1; k < n; ++k)
        {
            double term = coeffs_[k];

            // Product of (x - x_i) for i < k, excluding one factor
            for (size_t i = 0; i < k; ++i)
            {
                double product = 1.0;
                for (size_t j = 0; j < k; ++j)
                {
                    if (j != i)
                    {
                        product *= (x - x_[j]);
                    }
                }
                term *= product;
            }

            result += term;
        }

        return result;
    }
};

/**
 * @brief Convenience function for polynomial interpolation
 */
inline std::vector<double> interp1_polynomial(std::span<const double> x,
                                              std::span<const double> y,
                                              std::span<const double> x_new)
{
    return PolynomialInterpolator(x, y)(x_new);
}

} // namespace msl::interp

#endif // MSL_POLYNOMIAL_INTERPOLATOR_HPP