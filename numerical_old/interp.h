/**
**            qREST - Quick Response Evaluation for Safety Tagging
**     Institute of Engineering Mechanics, China Earthquake Administration
**
**                 Copyright 2024 - 2024 QLab, Dong Feiyue
**                          All Rights Reserved.
**
** Project: qREST
** File: \src\numerical_algorithm\interp.h
** -----
** File Created: Sunday, 14th July 2024 21:23:39
** Author: Dong Feiyue (donfeiyue@outlook.com)
** -----
** Last Modified: Saturday, 10th August 2024 22:29:07
** Modified By: Dong Feiyue (donfeiyue@outlook.com)
*/

// Description:
// 数值算法：插值算法

#ifndef NUMERICAL_ALGORITHM_INTERP_H
#define NUMERICAL_ALGORITHM_INTERP_H

// stdc++ headers
#include <vector>

// third-party library headers
#include "eigen3/Eigen/Core"
#include "gsl/gsl_interp.h"


namespace numerical_algorithm
{

// 插值算法类型
enum class InterpType
{
    Linear,
    CubicSpline,
    Akima,
    Steffen,
    Polynomial
};

// 插值算法类
class Interp
{
public:
    Interp() = default;

    // 由两个std::vector参数构造
    // @param x 插值点的x坐标
    // @param y 插值点的y坐标
    // @param interp_type 插值算法类型
    Interp(const std::vector<double> &x,
           const std::vector<double> &y,
           InterpType interp_type);

    // 由Eigen::MatrixXd参数构造
    // @param x 插值点的x坐标
    // @param y 插值点的y坐标
    // @param interp_type 插值算法类型
    Interp(const Eigen::VectorXd &x,
           const Eigen::VectorXd &y,
           InterpType interp_type);

    // 由插值方法参数构造
    // @param interp_type 插值算法类型
    explicit Interp(InterpType interp_type);

    // 析构函数
    ~Interp() = default;

    // 设置插值算法的x输入
    // @param x 插值点的x坐标
    void set_x(const std::vector<double> &x);

    // 设置插值算法的y输入
    // @param y 插值点的y坐标
    void set_y(const std::vector<double> &y);

    // 设置插值方法
    // @param interp_type 插值算法类型
    void set_interp_type(InterpType interp_type);

    // 计算单个点的插值结果
    // @param x_interp 待插值点的x坐标
    // @return 待插值点的插值结果y坐标
    double Interpolation(double x_interp) const;

    // 计算多个点的插值结果
    // @param x_interp 待插值点的x坐标
    // @return 待插值点的插值结果y坐标
    std::vector<double>
    Interpolation(const std::vector<double> &x_interp) const;

    // 计算矩阵的插值结果，按行插值
    // @param x 已知点的x坐标向量
    // @param y 已知点的y坐标矩阵
    // @param x_interp 待插值点的x坐标向量
    // @return 待插值点的插值结果y坐标矩阵
    std::vector<std::vector<double>>
    Interpolation(const std::vector<double> &x,
                  const std::vector<std::vector<double>> &y,
                  const std::vector<double> &x_interp);

private:
    // 插值算法的输入点
    std::vector<double> x_{}, y_{};
    // 插值算法类型
    const gsl_interp_type *interp_type_ = nullptr;
};

} // namespace numerical_algorithm

#endif // NUMERICAL_ALGORITHM_INTERP_H