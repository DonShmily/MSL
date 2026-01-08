/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: filter_apply.hpp (Corrected)
** -----
** 改进点：
** 1. 使用矩阵求解计算初始条件（与MATLAB一致）
** 2. 修正边界填充索引
** 3. 改进滤波器状态更新逻辑
*/

#ifndef MSL_FILTER_APPLY_HPP
#define MSL_FILTER_APPLY_HPP

#include <algorithm>
#include <cmath>
#include <span>
#include <stdexcept>
#include <vector>


namespace msl::signal
{

struct FilterCoefficients
{
    std::vector<double> b; // Numerator coefficients
    std::vector<double> a; // Denominator coefficients (a[0] = 1.0)

    FilterCoefficients() = default;
    FilterCoefficients(std::vector<double> num, std::vector<double> den)
        : b(std::move(num)), a(std::move(den))
    {}

    [[nodiscard]] size_t order() const
    {
        return std::max(a.size(), b.size()) - 1;
    }
};

// ============================================================================
// 简易矩阵求解（用于计算初始条件）
// ============================================================================

/**
 * @brief 简单的LU分解求解线性方程组 Ax = b
 */
inline std::vector<double>
solve_linear_system(const std::vector<std::vector<double>> &A,
                    const std::vector<double> &b)
{
    size_t n = b.size();
    if (n == 0)
        return {};

    // 创建增广矩阵的副本
    std::vector<std::vector<double>> M(n, std::vector<double>(n));
    std::vector<double> rhs(b);

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            M[i][j] = A[i][j];
        }
    }

    // 高斯消元法
    for (size_t k = 0; k < n; ++k)
    {
        // 寻找主元
        size_t max_row = k;
        double max_val = std::abs(M[k][k]);
        for (size_t i = k + 1; i < n; ++i)
        {
            if (std::abs(M[i][k]) > max_val)
            {
                max_val = std::abs(M[i][k]);
                max_row = i;
            }
        }

        // 交换行
        if (max_row != k)
        {
            std::swap(M[k], M[max_row]);
            std::swap(rhs[k], rhs[max_row]);
        }

        // 消元
        for (size_t i = k + 1; i < n; ++i)
        {
            double factor = M[i][k] / M[k][k];
            for (size_t j = k; j < n; ++j)
            {
                M[i][j] -= factor * M[k][j];
            }
            rhs[i] -= factor * rhs[k];
        }
    }

    // 回代求解
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i)
    {
        x[i] = rhs[i];
        for (size_t j = i + 1; j < n; ++j)
        {
            x[i] -= M[i][j] * x[j];
        }
        x[i] /= M[i][i];
    }

    return x;
}

// ============================================================================
// Forward-backward filtering (zero-phase)
// ============================================================================

/**
 * @brief 计算filtfilt的初始条件（改进版本，与MATLAB一致）
 *
 * 使用矩阵方程求解：sp * zi = b[1:] - b[0] * a[1:]
 * 其中sp是companion矩阵的变体
 */
inline std::vector<double> compute_filtfilt_zi(const FilterCoefficients &coeffs)
{
    size_t nfilt = std::max(coeffs.a.size(), coeffs.b.size());
    size_t n = nfilt - 1;

    if (n == 0)
    {
        return {}; // No initial conditions needed for order 0
    }

    // 构建稀疏矩阵 sp (与MATLAB filtfilt一致)
    // rows = [1:nfilt-1, 2:nfilt-1, 1:nfilt-2]
    // cols = [ones(1,nfilt-1), 2:nfilt-1, 2:nfilt-1]
    // data = [1+a(2), a(3:nfilt), ones(1,nfilt-2), -ones(1,nfilt-2)]

    std::vector<std::vector<double>> sp(n, std::vector<double>(n, 0.0));

    // 第一列: [1+a[1], a[2], ..., a[n]]
    sp[0][0] = 1.0 + coeffs.a[1];
    for (size_t i = 1; i < n; ++i)
    {
        sp[i][0] = (i + 1 < coeffs.a.size()) ? coeffs.a[i + 1] : 0.0;
    }

    // 次对角线: ones (从第2列到第n列)
    for (size_t i = 0; i < n - 1; ++i)
    {
        sp[i][i + 1] = 1.0;
    }

    // 下三角部分: -ones (从第2列到第n列)
    for (size_t i = 1; i < n; ++i)
    {
        sp[i][i] = -1.0;
    }

    // 右侧向量: b[1:] - b[0] * a[1:]
    std::vector<double> rhs(n);
    for (size_t i = 0; i < n; ++i)
    {
        double b_val = (i + 1 < coeffs.b.size()) ? coeffs.b[i + 1] : 0.0;
        double a_val = (i + 1 < coeffs.a.size()) ? coeffs.a[i + 1] : 0.0;
        rhs[i] = b_val - coeffs.b[0] * a_val;
    }

    // 求解线性系统
    return solve_linear_system(sp, rhs);
}

/**
 * @brief 使用初始条件应用滤波器
 */
inline std::vector<double> filter_with_zi(std::span<const double> signal,
                                          const FilterCoefficients &coeffs,
                                          std::vector<double> zi)
{
    if (coeffs.a.empty() || coeffs.b.empty())
    {
        throw std::invalid_argument("Filter coefficients cannot be empty");
    }

    // 归一化系数（确保a[0] = 1.0）
    std::vector<double> a = coeffs.a;
    std::vector<double> b = coeffs.b;

    if (std::abs(a[0] - 1.0) > 1e-10)
    {
        double a0 = a[0];
        for (auto &val : a)
            val /= a0;
        for (auto &val : b)
            val /= a0;
    }

    const size_t n = signal.size();
    size_t filter_order = std::max(a.size(), b.size());

    // 扩展系数向量
    a.resize(filter_order, 0.0);
    b.resize(filter_order, 0.0);
    zi.resize(filter_order, 0.0);

    std::vector<double> output(n);

    // Direct Form II 滤波实现（与版本2一致）
    for (size_t i = 0; i < n; ++i)
    {
        // 更新状态（从后向前）
        for (size_t order = filter_order - 1; order > 0; --order)
        {
            if (i >= order)
            {
                zi[order - 1] = b[order] * signal[i - order]
                                - a[order] * output[i - order] + zi[order];
            }
        }

        // 计算输出
        output[i] = b[0] * signal[i] + zi[0];
    }

    return output;
}

/**
 * @brief 零相位滤波（改进版本）
 *
 * 修正点：
 * 1. 使用正确的矩阵求解计算初始条件
 * 2. 修正边界填充的索引（与MATLAB一致）
 * 3. 改进滤波器实现
 *
 * @param signal Input signal
 * @param coeffs Filter coefficients
 * @return Zero-phase filtered signal
 */
inline std::vector<double> filtfilt(std::span<const double> signal,
                                    const FilterCoefficients &coeffs)
{
    const size_t len = signal.size();
    const size_t nfilt = std::max(coeffs.a.size(), coeffs.b.size());
    const size_t nfact = 3 * (nfilt - 1); // Edge transient length

    if (len <= nfact)
    {
        throw std::invalid_argument(
            "Input data too short! Must have length > 3 * filter_order");
    }

    // 计算初始条件（使用改进的矩阵求解方法）
    auto zi_base = compute_filtfilt_zi(coeffs);

    // 创建填充信号（修正索引以匹配版本2）
    std::vector<double> padded;
    padded.reserve(len + 2 * nfact);

    // 左填充: 2*signal[0] - signal[nfact:1:-1]
    // 对应版本2的 SubvectorReverse(input_signal, nfact, 1)
    for (size_t i = 0; i < nfact; ++i)
    {
        padded.push_back(2.0 * signal[0] - signal[nfact - i]);
    }

    // 原始信号
    padded.insert(padded.end(), signal.begin(), signal.end());

    // 右填充: 2*signal[end] - signal[end-2:end-nfact-1:-1]
    // 对应版本2的 SubvectorReverse(input_signal, len-2, len-nfact-1)
    for (size_t i = 0; i < nfact; ++i)
    {
        padded.push_back(2.0 * signal[len - 1] - signal[len - 2 - i]);
    }

    // 前向滤波
    std::vector<double> zi_forward = zi_base;
    for (auto &z : zi_forward)
    {
        z *= padded[0];
    }
    auto forward = filter_with_zi(padded, coeffs, zi_forward);

    // 反转
    std::reverse(forward.begin(), forward.end());

    // 后向滤波
    std::vector<double> zi_backward = zi_base;
    for (auto &z : zi_backward)
    {
        z *= forward[0];
    }
    auto backward = filter_with_zi(forward, coeffs, zi_backward);

    // 反转回原顺序
    std::reverse(backward.begin(), backward.end());

    // 提取有效部分
    std::vector<double> result(len);
    std::copy(backward.begin() + nfact,
              backward.begin() + nfact + len,
              result.begin());

    return result;
}

/**
 * @brief 零相位滤波（vector重载）
 */
inline std::vector<double> filtfilt(const std::vector<double> &signal,
                                    const FilterCoefficients &coeffs)
{
    return filtfilt(std::span<const double>(signal), coeffs);
}

} // namespace msl::signal

#endif // MSL_FILTER_APPLY_HPP