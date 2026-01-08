/**
**            qREST - Quick Response Evaluation for Safety Tagging
**     Institute of Engineering Mechanics, China Earthquake Administration
**
**                 Copyright 2024 - 2024 QLab, Dong Feiyue
**                          All Rights Reserved.
**
** Project: qREST
** File: \src\numerical_algorithm\math_function.h
** -----
** File Created: Thursday, 5th September 2024 00:06:50
** Author: Dong Feiyue (donfeiyue@outlook.com)
** -----
** Last Modified: Thursday, 5th September 2024 00:07:19
** Modified By: Dong Feiyue (donfeiyue@outlook.com)
*/

// Description: 定义了一些常用的数学函数，用于数值计算

#ifndef NUMERICAL_ALGORITHM_MATH_FUNCTION_H
#define NUMERICAL_ALGORITHM_MATH_FUNCTION_H

#define _USE_MATH_DEFINES

// stdc++ headers
#include <cmath>

namespace numerical_algorithm
{
// 计算同数量级的向上取整值
// @param x 待计算的数值
// @return 同数量级的向上取整值
// @note 例如：ceil_magnitude(0.0032) = 0.004，ceil_magnitude(163) = 200
inline double ceil_magnitude(double x)
{
    int sign = x < 0 ? x = std::abs(x), -1 : 1;
    double magnitude = std::pow(10, std::floor(std::log10(x)));
    return sign * std::ceil(x / magnitude) * magnitude;
}

// 计算同数量级的向下取整值
// @param x 待计算的数值
// @return 同数量级的向下取整值
// @note 例如：floor_magnitude(0.0032) = 0.003，floor_magnitude(163) = 100
inline double floor_magnitude(double x)
{
    int sign = x < 0 ? x = std::abs(x), -1 : 1;
    double magnitude = std::pow(10, std::floor(std::log10(x)));
    return sign * std::floor(x / magnitude) * magnitude;
}

} // namespace numerical_algorithm

#endif // NUMERICAL_ALGORITHM_MATH_FUNCTION_H