/**
**            qREST - Quick Response Evaluation for Safety Tagging
**     Institute of Engineering Mechanics, China Earthquake Administration
**
**                 Copyright 2024 - 2024 QLab, Dong Feiyue
**                          All Rights Reserved.
**
** Project: qREST
** File: \src\numerical_algorithm\filter.cpp
** -----
** File Created: Saturday, 13th July 2024 00:10:59
** Author: Dong Feiyue (donfeiyue@outlook.com)
** -----
** Last Modified: Saturday, 10th August 2024 22:28:26
** Modified By: Dong Feiyue (donfeiyue@outlook.com)
*/

// Description:
// 数值算法：filter滤波方法类的实现

// associated header
#include "filter.h"

// stdc++ headers
#include <stdexcept>
#include <vector>


namespace numerical_algorithm
{

// 矩阵滤波算法入口
std::vector<std::vector<double>>
Filter::Filtering(const std::vector<std::vector<double>> &input_signal)
{
    auto output_signal = std::vector<std::vector<double>>(
        input_signal.size(),
        std::vector<double>(input_signal.front().size(), 0.0));
    if (coefficients_a_.empty() || coefficients_b_.empty())
    {
        throw std::runtime_error("Filter coefficients are empty.");
    }
    if (input_signal.size() == 0 || input_signal.front().size() == 0)
    {
        throw std::runtime_error("Input signal is empty.");
    }
    if (coefficients_a_.front() == 0)
    {
        throw std::runtime_error("Filter coefficients are invalid.");
    }

    if (coefficients_a_.size() == 1)
    {
        fir_filter(input_signal, output_signal);
    }
    else
    {
        iir_filter(input_signal, output_signal);
    }
    return output_signal;
}

// 单列数据滤波算法入口
std::vector<double> Filter::Filtering(const std::vector<double> &input_signal)
{
    auto output_signal = std::vector<double>(input_signal.size(), 0.0);
    if (coefficients_a_.empty() || coefficients_b_.empty())
    {
        throw std::runtime_error("Filter coefficients are empty.");
    }
    if (input_signal.size() == 0)
    {
        throw std::runtime_error("Input signal is empty.");
    }
    if (coefficients_a_.front() == 0)
    {
        throw std::runtime_error("Filter coefficients are invalid.");
    }


    if (coefficients_a_.size() == 1)
    {
        fir_filter(input_signal, output_signal);
    }
    else
    {
        iir_filter(input_signal, output_signal);
    }
    return output_signal;
}

// FIR滤波算法
void Filter::fir_filter(const std::vector<std::vector<double>> &input_signal,
                        std::vector<std::vector<double>> &output_signal)
{
    output_signal = std::vector<std::vector<double>>(
        input_signal.size(),
        std::vector<double>(input_signal.front().size(), 0.0));
    for (std::size_t i = 0; i < input_signal.size(); ++i)
    {
        for (std::size_t j = 0; j < input_signal.front().size(); ++j)
        {
            for (std::size_t k = 0; k < coefficients_b_.size(); ++k)
            {
                if (j >= k)
                {
                    output_signal[i][j] +=
                        coefficients_b_[k] * input_signal[i][j - k];
                }
            }
        }
    }
}

// FIR滤波算法
void Filter::fir_filter(const std::vector<double> &input_signal,
                        std::vector<double> &output_signal)
{
    output_signal = std::vector<double>(input_signal.size(), 0.0);
    for (std::size_t j = 0; j < input_signal.size(); ++j)
    {
        for (std::size_t k = 0; k < coefficients_b_.size(); ++k)
        {
            if (j >= k)
            {
                output_signal[j] += coefficients_b_[k] * input_signal[j - k];
            }
        }
    }
}

// IIR滤波算法
void Filter::iir_filter(const std::vector<std::vector<double>> &input_signal,
                        std::vector<std::vector<double>> &output_signal)
{
    output_signal = std::vector<std::vector<double>>(
        input_signal.size(),
        std::vector<double>(input_signal.front().size(), 0.0));
    for (std::size_t i = 0; i < input_signal.size(); ++i)
    {
        for (std::size_t j = 0; j < input_signal.front().size(); ++j)
        {
            for (std::size_t k = 0; k < coefficients_b_.size(); ++k)
            {
                if (j >= k)
                {
                    output_signal[i][j] +=
                        coefficients_b_[k] * input_signal[i][j - k];
                }
            }
            for (std::size_t k = 1; k < coefficients_a_.size(); ++k)
            {
                if (j >= k)
                {
                    output_signal[i][j] -=
                        coefficients_a_[k] * output_signal[i][j - k];
                }
            }
        }
    }
}

// IIR滤波算法
void Filter::iir_filter(const std::vector<double> &input_signal,
                        std::vector<double> &output_signal)
{
    output_signal = std::vector<double>(input_signal.size(), 0.0);
    for (std::size_t j = 0; j < input_signal.size(); ++j)
    {
        for (std::size_t k = 0; k < coefficients_b_.size(); ++k)
        {
            if (j >= k)
            {
                output_signal[j] += coefficients_b_[k] * input_signal[j - k];
            }
        }
        for (std::size_t k = 1; k < coefficients_a_.size(); ++k)
        {
            if (j >= k)
            {
                output_signal[j] -= coefficients_a_[k] * output_signal[j - k];
            }
        }
    }
}

} // namespace numerical_algorithm