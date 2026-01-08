/**
**            qREST - Quick Response Evaluation for Safety Tagging
**     Institute of Engineering Mechanics, China Earthquake Administration
**
**                 Copyright 2024 - 2024 QLab, Dong Feiyue
**                          All Rights Reserved.
**
** Project: qREST
** File: \src\numerical_algorithm\integral.h
** -----
** File Created: Sunday, 14th July 2024 23:34:22
** Author: Dong Feiyue (donfeiyue@outlook.com)
** -----
** Last Modified: Saturday, 10th August 2024 22:28:52
** Modified By: Dong Feiyue (donfeiyue@outlook.com)
*/

// Description:
// 数值算法：积分算法

#ifndef NUMERICAL_ALGORITHM_INTEGRAL_H
#define NUMERICAL_ALGORITHM_INTEGRAL_H

// stdc++ headers
#include <vector>

namespace numerical_algorithm
{
// std::vector<double>梯形积分算法
// @param input 输入数据
// @param output 积分结果
// @param dx 积分步长
inline std::vector<double> Cumtrapz(const std::vector<double> &input, double dx)
{
    std::vector<double> output(input.size(), 0.0);
    for (size_t i = 1; i < input.size(); ++i)
    {
        output[i] = output[i - 1] + 0.5 * (input[i] + input[i - 1]) * dx;
    }
    return output;
}

// std::vector<std::vector<double>> 梯形积分算法（按列积分）
// @param input 输入数据矩阵
// @param output 积分结果矩阵
// @param dx 积分步长
inline std::vector<std::vector<double>>
Cumtrapz(const std::vector<std::vector<double>> &input, double dx)
{
    std::vector<std::vector<double>> output(
        input.size(), std::vector<double>(input[0].size(), 0.0));
    for (size_t i = 0; i < input.size(); ++i)
    {
        for (size_t j = 1; j < input[i].size(); ++j)
        {
            output[i][j] =
                output[i][j - 1] + 0.5 * (input[i][j] + input[i][j - 1]) * dx;
        }
    }
    return output;
}

} // namespace numerical_algorithm

#endif // NUMERERICAL_ALGORITHM_INTEGRAL_H