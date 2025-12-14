/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: polynominal.hpp
** -----
** File Created: Saturday, 13th December 2025 22:56:46
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:05:05
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_POLYNOMIAL_HPP
#define MSL_POLYNOMIAL_HPP

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

#include "matrix/martrix_decompose.hpp"
#include "matrix/matrix_operation.hpp"
#include "matrix/real_matrix_owned.hpp"

namespace msl::polynomial
{
// Polynomial evaluation
class Polynomial
{
private:
    std::vector<double> coeffs_{};

public:
    Polynomial() = default;
    Polynomial(const std::vector<double> &coeffs) : coeffs_(coeffs) {}
    Polynomial(std::span<const double> y, std::size_t n = 0)
    {
        std::vector<double> x(y.size());
        std::iota(x.begin(), x.end(), 0.0);
        calc_coefficients(x, y, n);
    }
    Polynomial(std::span<const double> x,
               std::span<const double> y,
               std::size_t n = 0)
    {
        calc_coefficients(x, y, n);
    }

    double operator()(double x) const { return evaluate(x); }

    std::vector<double>
    operator()(const std::span<const double> &x_values) const
    {
        std::vector<double> results(x_values.size());
        for (std::size_t i = 0; i < x_values.size(); ++i)
        {
            results[i] = evaluate(x_values[i]);
        }
        return results;
    }

    std::vector<double> coefficients() const { return coeffs_; }

private:
    double evaluate(double x) const
    {
        if (coeffs_.empty())
        {
            throw std::runtime_error("Polynomial coefficients are not set");
        }

        double result = 0.0;
        double x_pow = 1.0; // x^0

        for (const auto &c : coeffs_)
        {
            result += c * x_pow;
            x_pow *= x;
        }

        return result;
    }

    void calc_coefficients(std::span<const double> x,
                           std::span<const double> y,
                           std::size_t n = 0)
    {
        if (x.size() != y.size())
        {
            throw std::invalid_argument("x and y must have same size");
        }

        if (x.size() < n + 1)
        {
            throw std::invalid_argument(
                "Not enough points to fit the polynomial");
        }

        if (n == 0)
        {
            coeffs_.resize(1);
            coeffs_[0] = std::accumulate(y.begin(), y.end(), 0.0)
                         / static_cast<double>(y.size());
            return;
        }

        // least squares fitting
        matrix::matrixd A(x.size(), n + 1);
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            double x_pow = 1.0;
            for (std::size_t j = 0; j <= n; ++j)
            {
                A(i, j) = x_pow;
                x_pow *= x[i];
            }
        }

        auto qr_result = matrix::qr(A);
        auto &Q = qr_result[0];
        auto &R = qr_result[1];

        auto Qt_y =
            matrix::transpose(Q) * matrix::real_matrix_owned(y.size(), 1, y);
        Qt_y = Qt_y.submatrix(0, n + 1, 0, 1);
        matrix::real_matrix_owned coeffs_mat =
            matrix::real_matrix_owned(n + 1, 1);
        matrix::matrixd R_upper(n + 1, n + 1, 0.0);
        for (std::size_t i = 0; i <= n; ++i)
        {
            for (std::size_t j = i; j <= n; ++j)
            {
                R_upper(i, j) = R(i, j);
            }
        }

        auto c = matrix::inverse(R_upper) * Qt_y;

        coeffs_.resize(n + 1);
        std::copy(c.data(), c.data() + c.size(), coeffs_.begin());
    }
};
}; // namespace msl::polynomial

#endif // MSL_POLYNOMIAL_HPP
