/**
**            qREST - Quick Response Evaluation for Safety Tagging
**     Institute of Engineering Mechanics, China Earthquake Administration
**
**                 Copyright 2024 - 2024 QLab, Dong Feiyue
**                          All Rights Reserved.
**
** Project: qREST
** File: \src\numerical_algorithm\basic_filtering.h
** -----
** File Created: Friday, 12th July 2024 23:32:48
** Author: Dong Feiyue (donfeiyue@outlook.com)
** -----
** Last Modified: Saturday, 10th August 2024 22:28:06
** Modified By: Dong Feiyue (donfeiyue@outlook.com)
*/

// Description:
// 数值算法：滤波算法基类

#ifndef NUMERICAL_ALGORITHM_BASIC_FILTERING_H
#define NUMERICAL_ALGORITHM_BASIC_FILTERING_H

// stdc++ headers
#include <vector>


namespace numerical_algorithm
{

// 滤波函数类型
enum class FilterFunction
{
    filter,
    filtfilt
};

// 滤波算法基类
class BasicFiltering
{
public:
    // 默认构造函数
    BasicFiltering() = default;

    // 析构函数
    virtual ~BasicFiltering() = default;

    // 矩阵滤波算法入口
    // @param input_signal 输入信号矩阵
    // @return 滤波后的信号矩阵
    virtual std::vector<std::vector<double>>
    Filtering(const std::vector<std::vector<double>> &input_signa) = 0;

    // 单列数据滤波算法入口
    // @param input_signal 输入信号
    // @return 滤波后的信号
    virtual std::vector<double>
    Filtering(const std::vector<double> &input_signal) = 0;

protected:
};

} // namespace numerical_algorithm

#endif // NUMERICAL_ALGORITHM_BASIC_FILTERING_H
