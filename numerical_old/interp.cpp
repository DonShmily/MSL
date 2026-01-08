/**
**            qREST - Quick Response Evaluation for Safety Tagging
**     Institute of Engineering Mechanics, China Earthquake Administration
**
**                 Copyright 2024 - 2024 QLab, Dong Feiyue
**                          All Rights Reserved.
**
** Project: qREST
** File: \src\numerical_algorithm\interp.cpp
** -----
** File Created: Sunday, 14th July 2024 21:23:47
** Author: Dong Feiyue (donfeiyue@outlook.com)
** -----
** Last Modified: Saturday, 10th August 2024 22:28:59
** Modified By: Dong Feiyue (donfeiyue@outlook.com)
*/

// Description:
// 数值算法：插值算法

// associated header
#include "interp.h"

// stdc++ headers
#include <stdexcept>
#include <vector>

// third-party library headers
#include <gsl/gsl_interp.h>
#include "eigen3/Eigen/Dense"
#include "gsl/gsl_spline.h"


namespace numerical_algorithm
{
// 由两个std::vector参数构造
Interp::Interp(const std::vector<double> &x,
               const std::vector<double> &y,
               InterpType interp_type)
    : x_(x), y_(y)
{
    if (x.size() != y.size())
    {
        throw std::invalid_argument("The size of x and y must be equal.");
    }
    switch (interp_type)
    {
        case InterpType::Linear:
            interp_type_ = gsl_interp_linear;
            break;
        case InterpType::CubicSpline:
            interp_type_ = gsl_interp_cspline;
            break;
        case InterpType::Akima:
            interp_type_ = gsl_interp_akima;
            break;
        case InterpType::Steffen:
            interp_type_ = gsl_interp_steffen;
            break;
        case InterpType::Polynomial:
            interp_type_ = gsl_interp_polynomial;
            break;
        default:
            interp_type_ = gsl_interp_linear;
            break;
    }
}

// 由Eigen::MatrixXd参数构造
Interp::Interp(const Eigen::VectorXd &x,
               const Eigen::VectorXd &y,
               InterpType interp_type)
{
    if (x.size() != y.size())
    {
        throw std::invalid_argument("The size of x and y must be equal.");
    }
    x_.resize(x.size());
    y_.resize(y.size());
    for (int i = 0; i < x.size(); ++i)
    {
        x_[i] = x(i);
        y_[i] = y(i);
    }
    switch (interp_type)
    {
        case InterpType::Linear:
            interp_type_ = gsl_interp_linear;
            break;
        case InterpType::CubicSpline:
            interp_type_ = gsl_interp_cspline;
            break;
        case InterpType::Akima:
            interp_type_ = gsl_interp_akima;
            break;
        case InterpType::Steffen:
            interp_type_ = gsl_interp_steffen;
            break;
        case InterpType::Polynomial:
            interp_type_ = gsl_interp_polynomial;
            break;
        default:
            interp_type_ = gsl_interp_linear;
            break;
    }
}

// 由插值方法参数构造
Interp::Interp(InterpType interp_type)
{
    switch (interp_type)
    {
        case InterpType::Linear:
            interp_type_ = gsl_interp_linear;
            break;
        case InterpType::CubicSpline:
            interp_type_ = gsl_interp_cspline;
            break;
        case InterpType::Akima:
            interp_type_ = gsl_interp_akima;
            break;
        case InterpType::Steffen:
            interp_type_ = gsl_interp_steffen;
            break;
        case InterpType::Polynomial:
            interp_type_ = gsl_interp_polynomial;
            break;
        default:
            interp_type_ = gsl_interp_linear;
            break;
    }
}

// 设置插值算法的x输入
void Interp::set_x(const std::vector<double> &x) { x_ = x; }

// 设置插值算法的y输入
void Interp::set_y(const std::vector<double> &y) { y_ = y; }

// 设置插值方法
void Interp::set_interp_type(InterpType interp_type)
{
    switch (interp_type)
    {
        case InterpType::Linear:
            interp_type_ = gsl_interp_linear;
            break;
        case InterpType::CubicSpline:
            interp_type_ = gsl_interp_cspline;
            break;
        case InterpType::Akima:
            interp_type_ = gsl_interp_akima;
            break;
        case InterpType::Steffen:
            interp_type_ = gsl_interp_steffen;
            break;
        case InterpType::Polynomial:
            interp_type_ = gsl_interp_polynomial;
            break;
        default:
            interp_type_ = gsl_interp_linear;
            break;
    }
}

// 计算单个点的插值结果
double Interp::Interpolation(double x_interp) const
{
    if (x_.empty() || y_.empty())
    {
        throw std::domain_error("Input point is not set.");
    }
    if (x_interp < x_.front())
    {
        return (y_.front() - y_[1]) / (x_.front() - x_[1])
                   * (x_interp - x_.front())
               + y_.front();
    }
    else if (x_interp > x_.back())
    {
        return (y_.back() - y_[y_.size() - 2]) / (x_.back() - x_[x_.size() - 2])
                   * (x_interp - x_.back())
               + y_.back();
    }
    gsl_interp_accel *accel = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(interp_type_, x_.size());
    gsl_spline_init(spline, x_.data(), y_.data(), x_.size());
    double y_interp = gsl_spline_eval(spline, x_interp, accel);
    gsl_spline_free(spline);
    gsl_interp_accel_free(accel);
    return y_interp;
}

// 计算多个点的插值结果
std::vector<double>
Interp::Interpolation(const std::vector<double> &x_interp) const
{
    std::vector<double> y_interp(x_interp.size());
    if (x_.empty() || y_.empty())
    {
        throw std::domain_error("Input point is not set.");
    }
    gsl_interp_accel *accel = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(interp_type_, x_.size());
    gsl_spline_init(spline, x_.data(), y_.data(), x_.size());
    for (int i = 0; i < x_interp.size(); ++i)
    {
        if (x_interp[i] < x_.front())
        {
            y_interp[i] = (y_.front() - y_[1]) / (x_.front() - x_[1])
                              * (x_interp[i] - x_.front())
                          + y_.front();
        }
        else if (x_interp[i] > x_.back())
        {
            y_interp[i] = (y_.back() - y_[y_.size() - 2])
                              / (x_.back() - x_[x_.size() - 2])
                              * (x_interp[i] - x_.back())
                          + y_.back();
        }
        else
        {
            y_interp[i] = gsl_spline_eval(spline, x_interp[i], accel);
        }
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(accel);
    return y_interp;
}

// 计算矩阵的插值结果
std::vector<std::vector<double>>
Interp::Interpolation(const std::vector<double> &x,
                      const std::vector<std::vector<double>> &y,
                      const std::vector<double> &x_interp)
{
    if (x.empty() || y.empty())
    {
        throw std::domain_error("Input point is empty.");
    }
    std::vector<std::vector<double>> y_interp(
        x_interp.size(), std::vector<double>(y.front().size()));
    if (x.size() != y.size())
    {
        throw std::invalid_argument("The size of x and y must be equal.");
    }
    x_ = x;
    y_.resize(x_.size());
    gsl_interp_accel *accel = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(interp_type_, x_.size());
    for (int i = 0; i < y.front().size(); ++i)
    {
        for (int j = 0; j < x_.size(); ++j)
        {
            y_[j] = y[j][i];
        }
        gsl_spline_init(spline, x_.data(), y_.data(), x_.size());
        for (int j = 0; j < x_interp.size(); ++j)
        {
            if (x_interp[j] < x_.front())
            {
                y_interp[j][i] = (y_.front() - y_[1]) / (x_.front() - x_[1])
                                     * (x_interp[j] - x_.front())
                                 + y_.front();
            }
            else if (x_interp[j] > x_.back())
            {
                y_interp[j][i] = (y_.back() - y_[y_.size() - 2])
                                     / (x_.back() - x_[x_.size() - 2])
                                     * (x_interp[j] - x_.back())
                                 + y_.back();
            }
            else
                y_interp[j][i] = gsl_spline_eval(spline, x_interp[j], accel);
        }
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(accel);
    return y_interp;
}

} // namespace numerical_algorithm