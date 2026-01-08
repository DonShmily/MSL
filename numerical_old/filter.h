/**
**            qREST - Quick Response Evaluation for Safety Tagging
**     Institute of Engineering Mechanics, China Earthquake Administration
**
**                 Copyright 2024 - 2024 QLab, Dong Feiyue
**                          All Rights Reserved.
**
** Project: qREST
** File: \src\numerical_algorithm\filter.h
** -----
** File Created: Saturday, 13th July 2024 00:10:50
** Author: Dong Feiyue (donfeiyue@outlook.com)
** -----
** Last Modified: Saturday, 10th August 2024 22:28:33
** Modified By: Dong Feiyue (donfeiyue@outlook.com)
*/

// Description:
// 数值算法：filter滤波方法类

#ifndef NUMERICAL_ALGORITHM_FILTER_H
#define NUMERICAL_ALGORITHM_FILTER_H

// stdc++ headers
#include <vector>

// project headers
#include "basic_filtering.h"
#include "butterworth_filter_design.h"


namespace numerical_algorithm
{

// Filter滤波方法类
class Filter : public BasicFiltering
{
public:
    // 默认构造函数
    Filter() = default;

    // 由两个std::vector参数构造
    // @param coefficients_a 滤波器分母系数
    // @param coefficients_b 滤波器分子系数
    Filter(std::vector<double> coefficients_a,
           std::vector<double> coefficients_b)
        : coefficients_a_(coefficients_a), coefficients_b_(coefficients_b)
    {}

    // 由一个std::vector作为参数b构造，a默认为1
    // @param coefficients_b 滤波器分子系数
    explicit Filter(std::vector<double> coefficients_b)
        : coefficients_a_(1, 1.0), coefficients_b_(coefficients_b)
    {}

    // 由butterworth滤波器设计构造
    // @param filter_design butterworth滤波器设计类
    explicit Filter(const ButterworthFilterDesign &filter_design)
    {
        filter_design.get_filter_coefficients(coefficients_a_, coefficients_b_);
    }

    // 析构函数
    ~Filter() = default;

    // 矩阵滤波算法的入口
    // @param input_signal 输入信号矩阵
    // @return 滤波后的信号矩阵
    std::vector<std::vector<double>>
    Filtering(const std::vector<std::vector<double>> &input_signal) override;

    // 单列数据滤波算法的入口
    // @param input_signal 输入信号
    // @return 滤波后的信号
    std::vector<double>
    Filtering(const std::vector<double> &input_signal) override;

private:
    // filter滤波方法参数
    std::vector<double> coefficients_a_, coefficients_b_;

    // 有限脉冲响应滤波器（矩阵）
    void fir_filter(const std::vector<std::vector<double>> &input_signal,
                    std::vector<std::vector<double>> &output_signal);

    // 有限脉冲响应滤波器（单列）
    void fir_filter(const std::vector<double> &input_signal,
                    std::vector<double> &output_signal);

    // 无限脉冲响应滤波器（矩阵）
    void iir_filter(const std::vector<std::vector<double>> &input_signal,
                    std::vector<std::vector<double>> &output_signal);

    // 无限脉冲响应滤波器（单列）
    void iir_filter(const std::vector<double> &input_signal,
                    std::vector<double> &output_signal);
};

} // namespace numerical_algorithm

#endif // NUMERICAL_ALGORITHM_FILTER_H