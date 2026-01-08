/**
**            qREST - Quick Response Evaluation for Safety Tagging
**     Institute of Engineering Mechanics, China Earthquake Administration
**
**                 Copyright 2024 - 2024 QLab, Dong Feiyue
**                          All Rights Reserved.
**
** Project: qREST
** File: \src\numerical_algorithm\filtfilt.h
** -----
** File Created: Saturday, 13th July 2024 23:37:23
** Author: Dong Feiyue (donfeiyue@outlook.com)
** -----
** Last Modified: Saturday, 10th August 2024 22:28:46
** Modified By: Dong Feiyue (donfeiyue@outlook.com)
*/

// Description：
// 数值算法：filtfilt滤波方法类

#ifndef NUMERICAL_ALGORITHM_FILTFILT_H
#define NUMERICAL_ALGORITHM_FILTFILT_H

// stdc++ headers
#include <algorithm>
#include <vector>

// third-party headers

// project headers
#include "basic_filtering.h"
#include "butterworth_filter_design.h"

namespace numerical_algorithm
{

// filtfilt滤波方法类
class FiltFilt : public BasicFiltering
{
public:
    // 默认构造函数
    FiltFilt() = default;

    // 由两个std::vector参数构造
    // @param coefficients_a 滤波器分母系数
    // @param coefficients_b 滤波器分子系数
    FiltFilt(std::vector<double> coefficients_a,
             std::vector<double> coefficients_b)
        : coefficients_a_(coefficients_a), coefficients_b_(coefficients_b)
    {}

    // 由一个std::vector作为参数b构造，a默认为1
    // @param coefficients_b 滤波器分子系数
    explicit FiltFilt(std::vector<double> coefficients_b)
        : coefficients_a_(1, 1.0), coefficients_b_(coefficients_b)
    {}

    // 由ButterworthFilterDesign参数构造
    // @param filter_design butterworth滤波器设计类
    FiltFilt(const ButterworthFilterDesign &filter_design)
    {
        filter_design.get_filter_coefficients(coefficients_a_, coefficients_b_);
    }

    // 析构函数
    ~FiltFilt() = default;

    // 设置滤波方法参数
    // @param coefficients_a 滤波器分母系数
    // @param coefficients_b 滤波器分子系数
    void set_coefficients(std::vector<double> coefficients_a,
                          std::vector<double> coefficients_b)
    {
        coefficients_a_ = coefficients_a;
        coefficients_b_ = coefficients_b;
    }

    // 设置滤波方法参数
    // @param coefficients_b 滤波器分子系数
    void set_coefficients(std::vector<double> coefficients_b)
    {
        coefficients_a_ = std::vector<double>(1, 1.0);
        coefficients_b_ = coefficients_b;
    }

    // 设置滤波方法参数
    // @param filter_design butterworth滤波器设计类
    void set_coefficients(const ButterworthFilterDesign &filter_design)
    {
        filter_design.get_filter_coefficients(coefficients_a_, coefficients_b_);
    }

    // 矩阵滤波算法入口
    // @param input_signal 输入信号矩阵
    // @return 滤波后的信号矩阵
    std::vector<std::vector<double>>
    Filtering(const std::vector<std::vector<double>> &input_signal) override;

    // 单列数据滤波算法入口
    // @param input_signal 输入信号
    // @return 滤波后的信号
    std::vector<double>
    Filtering(const std::vector<double> &input_signal) override;

private:
    // filtfilt滤波方法参数
    std::vector<double> coefficients_a_, coefficients_b_;

    // filtfilt滤波算法的filter函数
    void filter(const std::vector<double> &input_signal,
                std::vector<double> &output_signal,
                std::vector<double> zi);

    // filtfilt算法的功能函数
    inline void
    AddIndexRange(std::vector<int> &indices, int beg, int end, int inc = 1)
    {
        for (int i = beg; i <= end; i += inc)
        {
            indices.push_back(i);
        }
    }

    inline void
    AddIndexConst(std::vector<int> &indices, int value, size_t numel)
    {
        indices.resize(numel + indices.size(), value);
    }

    inline void AppendVector(std::vector<double> &vec,
                             const std::vector<double> &tail)
    {
        vec.insert(vec.end(), tail.begin(), tail.end());
    }

    inline std::vector<double>
    SubvectorReverse(const std::vector<double> &vec, int idx_end, int idx_start)
    {
        std::vector<double> result(&vec[idx_start], &vec[idx_end + 1]);
        std::reverse(result.begin(), result.end());
        return result;
    }

    inline int MaxVal(const std::vector<int> &vec)
    {
        return *std::max_element(vec.begin(), vec.end());
    }
};

} // namespace numerical_algorithm

#endif // NUMERICAL_ALGORITHM_FILTFILT_H