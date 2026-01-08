/**
**            qREST - Quick Response Evaluation for Safety Tagging
**     Institute of Engineering Mechanics, China Earthquake Administration
**
**                 Copyright 2024 - 2024 QLab, Dong Feiyue
**                          All Rights Reserved.
**
** Project: qREST
** File: \src\numerical_algorithm\vector_calculation.h
** -----
** File Created: Monday, 15th July 2024 17:00:30
** Author: Dong Feiyue (donfeiyue@outlook.com)
** -----
** Last Modified: Saturday, 10th August 2024 22:29:14
** Modified By: Dong Feiyue (donfeiyue@outlook.com)
*/

// Description:
// 数值算法：向量计算，拓展std::vector的操作

#ifndef NUMERICAL_ALGORITHM_VECTOR_CALCULATION_H
#define NUMERICAL_ALGORITHM_VECTOR_CALCULATION_H

// stdc++ headers
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <stdexcept>
#include <vector>

// third-party headers
#include "fftw3.h"

namespace numerical_algorithm
{

// 寻找向量中的极值点
// @param vec 向量
// @param threshold 阈值
// @return 极值点的索引
inline std::vector<std::size_t> FindPeaks(const std::vector<double> &vec,
                                          double threshold = 0.0)
{
    std::vector<std::size_t> peaks;
    for (std::size_t i = 1; i < vec.size() - 1; ++i)
    {
        if (vec[i] > vec[i - 1] && vec[i] > vec[i + 1] && vec[i] > threshold)
        {
            peaks.push_back(i);
        }
    }
    return peaks;
}

// 向量计算
// @param vector1 向量1
// @param vector2 向量2
// @param operation 操作符
// @return 计算结果
inline std::vector<double> VectorOperation(const std::vector<double> &vector1,
                                           const std::vector<double> &vector2,
                                           const char operation)
{
    auto vec_size = vector1.size();
    std::vector<double> result(vec_size);
    switch (operation)
    {
        case '+':
            for (std::size_t i = 0; i < vec_size; ++i)
            {
                result[i] = vector1[i] + vector2[i];
            }
            break;
        case '-':
            for (std::size_t i = 0; i < vec_size; ++i)
            {
                result[i] = vector1[i] - vector2[i];
            }
            break;
        case '*':
            for (std::size_t i = 0; i < vec_size; ++i)
            {
                result[i] = vector1[i] * vector2[i];
            }
            break;
        case '/':
            for (std::size_t i = 0; i < vec_size; ++i)
            {
                result[i] = vector1[i] / vector2[i];
            }
            break;
        default:
            throw std::invalid_argument("Invalid operation");
            break;
    }
    return result;
}

// 向量计算
// @param vector1 向量1
// @param scalar 标量
// @param operation 操作符
// @return 计算结果
inline std::vector<double> VectorOperation(const std::vector<double> &vector1,
                                           const double scalar,
                                           const char operation)
{
    auto vec_size = vector1.size();
    std::vector<double> result(vec_size);
    switch (operation)
    {
        case '+':
            for (std::size_t i = 0; i < vec_size; ++i)
            {
                result[i] = vector1[i] + scalar;
            }
            break;
        case '-':
            for (std::size_t i = 0; i < vec_size; ++i)
            {
                result[i] = vector1[i] - scalar;
            }
            break;
        case '*':
            for (std::size_t i = 0; i < vec_size; ++i)
            {
                result[i] = vector1[i] * scalar;
            }
            break;
        case '/':
            for (std::size_t i = 0; i < vec_size; ++i)
            {
                result[i] = vector1[i] / scalar;
            }
            break;
        default:
            throw std::invalid_argument("Invalid operation");
            break;
    }
    return result;
}

// 寻找向量绝对值最大的元素
// @param input_vector 输入向量
// @return 绝对值最大的元素
inline double FindMaxAbs(const std::vector<double> &input_vector)
{
    return std::abs(*std::max_element(
        input_vector.begin(), input_vector.end(), [](double a, double b) {
            return std::abs(a) < std::abs(b);
        }));
}

// 向量归一化
// @param input_vector 输入向量
// @return 归一化后的向量
inline std::vector<double> Normalize(const std::vector<double> &input_vector)
{
    auto max_abs = FindMaxAbs(input_vector);
    return VectorOperation(input_vector, max_abs, '/');
}

// 计算给定移位参数的互相关值
// @param vector_x 向量x
// @param vector_y 向量y
// @param shift 移位参数
inline double CrossCorrelationAtShift(const std::vector<double> &vector_x,
                                      const std::vector<double> &vector_y,
                                      int shift)
{
    int n_a = vector_x.size();
    int n_b = vector_y.size();
    double result = 0.0;
    for (int i = 0; i < n_a; ++i)
    {
        int j = i - shift;
        if (j >= 0 && j < n_b)
        {
            result += vector_x[i] * vector_y[j];
        }
    }
    return result;
}

// 计算给定移位参数的互相关值，并作归一化处理
// @param vector_x 向量x
// @param vector_y 向量y
// @param shift 移位参数
inline double
CrossCorrelationAtShiftNormalized(const std::vector<double> &vector_x,
                                  const std::vector<double> &vector_y,
                                  int shift)
{
    return CrossCorrelationAtShift(vector_x, vector_y, shift)
           / std::sqrt(CrossCorrelationAtShift(vector_x, vector_x, 0)
                       * CrossCorrelationAtShift(vector_y, vector_y, 0));
}

// 向量互相关序列
// @param vector_x 向量x
// @param vector_y 向量y
// @return 互相关序列
inline std::vector<double> CrossCorrelation(const std::vector<double> &vector_x,
                                            const std::vector<double> &vector_y)
{
    int n_a = vector_x.size();
    int n_b = vector_y.size();
    int n = n_a + n_b - 1;
    std::vector<double> result(n, 0.0);

    for (int shift = -(n_b - 1); shift <= n_a - 1; ++shift)
    {
        result[shift + (n_b - 1)] =
            CrossCorrelationAtShift(vector_x, vector_y, shift);
    }
    return result;
}

// 向量互相关序列，并作归一化处理
// @param vector_x 向量x
// @param vector_y 向量y
// @return 归一化后的互相关序列
inline std::vector<double>
CrossCorrelationNormalized(const std::vector<double> &vector_x,
                           const std::vector<double> &vector_y)
{
    auto result = CrossCorrelation(vector_x, vector_y);
    double norm = std::sqrt(CrossCorrelationAtShift(vector_x, vector_x, 0)
                            * CrossCorrelationAtShift(vector_y, vector_y, 0));
    return VectorOperation(result, norm, '/');
}

// 向量自相关序列
// @param vector_x 向量x
// @return 自相关序列
inline std::vector<double> AutoCorrelation(const std::vector<double> &vector_x)
{
    return CrossCorrelation(vector_x, vector_x);
}

// 向量卷积序列
// @param vector_x 向量x
// @param vector_y 向量y
// @return 卷积序列
inline std::vector<double> Convolution(const std::vector<double> &vector_x,
                                       const std::vector<double> &vector_y)
{
    int n_a = vector_x.size();
    int n_b = vector_y.size();
    int n = n_a + n_b - 1;
    std::vector<double> result(n, 0.0);

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n_b; ++j)
        {
            if ((i - j) >= 0 && (i - j) < n_a)
            {
                result[i] += vector_x[i - j] * vector_y[j];
            }
        }
    }
    return result;
}

// 向量Fourierr变换
// @param input_vector 输入向量
// @return Fourier变换后的复数向量
inline std::vector<std::complex<double>>
FourierTransform(const std::vector<double> &input_vector)
{
    int n = input_vector.size();
    std::vector<std::complex<double>> result(n, 0.0);
    fftw_plan plan =
        fftw_plan_dft_r2c_1d(n,
                             const_cast<double *>(input_vector.data()),
                             reinterpret_cast<fftw_complex *>(result.data()),
                             FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    return result;
}


} // namespace numerical_algorithm

#endif // NUMERICAL_ALGORITHM_VECTOR_CALCULATION_H