/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: filter_apply.hpp (Final Corrected Version)
** -----
** 完全按照版本2的逻辑重写，确保与MATLAB一致
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
// 矩阵运算辅助函数
// ============================================================================

inline std::vector<std::vector<double>>
matrix_inverse(const std::vector<std::vector<double>> &mat)
{
    size_t n = mat.size();
    std::vector<std::vector<double>> A = mat;
    std::vector<std::vector<double>> inv(n, std::vector<double>(n, 0.0));

    // 初始化为单位矩阵
    for (size_t i = 0; i < n; ++i)
    {
        inv[i][i] = 1.0;
    }

    // 高斯-约旦消元法
    for (size_t i = 0; i < n; ++i)
    {
        // 寻找主元
        size_t max_row = i;
        double max_val = std::abs(A[i][i]);
        for (size_t k = i + 1; k < n; ++k)
        {
            if (std::abs(A[k][i]) > max_val)
            {
                max_val = std::abs(A[k][i]);
                max_row = k;
            }
        }

        if (max_row != i)
        {
            std::swap(A[i], A[max_row]);
            std::swap(inv[i], inv[max_row]);
        }

        // 归一化当前行
        double pivot = A[i][i];
        for (size_t j = 0; j < n; ++j)
        {
            A[i][j] /= pivot;
            inv[i][j] /= pivot;
        }

        // 消元
        for (size_t k = 0; k < n; ++k)
        {
            if (k != i)
            {
                double factor = A[k][i];
                for (size_t j = 0; j < n; ++j)
                {
                    A[k][j] -= factor * A[i][j];
                    inv[k][j] -= factor * inv[i][j];
                }
            }
        }
    }

    return inv;
}

inline std::vector<double>
matrix_vector_multiply(const std::vector<std::vector<double>> &mat,
                       const std::vector<double> &vec)
{
    size_t rows = mat.size();
    size_t cols = vec.size();
    std::vector<double> result(rows, 0.0);

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            result[i] += mat[i][j] * vec[j];
        }
    }

    return result;
}

// ============================================================================
// Forward-backward filtering (zero-phase)
// ============================================================================

/**
 * @brief 计算filtfilt的初始条件（完全按照版本2的逻辑）
 */
inline std::vector<double> compute_filtfilt_zi(const FilterCoefficients &coeffs)
{
    size_t nfilt = std::max(coeffs.a.size(), coeffs.b.size());
    size_t n = nfilt - 1;

    if (n == 0)
    {
        return {};
    }

    // 扩展系数向量
    std::vector<double> a = coeffs.a;
    std::vector<double> b = coeffs.b;
    a.resize(nfilt, 0.0);
    b.resize(nfilt, 0.0);

    // 构建稀疏矩阵sp - 完全按照版本2的方式
    std::vector<int> rows, cols;
    std::vector<double> data;

    // rows = [0:n-1, 1:n-1, 0:n-2]
    for (size_t i = 0; i < n; ++i)
        rows.push_back(i);
    if (n > 1)
    {
        for (size_t i = 1; i < n; ++i)
            rows.push_back(i);
        for (size_t i = 0; i < n - 1; ++i)
            rows.push_back(i);
    }

    // cols = [0, 0, ..., 0, 1:n-1, 1:n-1]
    for (size_t i = 0; i < n; ++i)
        cols.push_back(0);
    if (n > 1)
    {
        for (size_t i = 1; i < n; ++i)
            cols.push_back(i);
        for (size_t i = 1; i < n; ++i)
            cols.push_back(i);
    }

    // data = [1+a[1], a[2:n], ones(1,n-1), -ones(1,n-1)]
    data.push_back(1.0 + a[1]);
    for (size_t i = 2; i < nfilt; ++i)
    {
        data.push_back(a[i]);
    }
    if (n > 1)
    {
        for (size_t i = 0; i < n - 1; ++i)
            data.push_back(1.0);
        for (size_t i = 0; i < n - 1; ++i)
            data.push_back(-1.0);
    }

    // 构建完整矩阵
    std::vector<std::vector<double>> sp(n, std::vector<double>(n, 0.0));
    for (size_t k = 0; k < rows.size(); ++k)
    {
        sp[rows[k]][cols[k]] += data[k]; // 注意：使用+=因为可能有重复位置
    }

    // 计算右侧向量: b[1:] - b[0] * a[1:]
    std::vector<double> rhs(n);
    for (size_t i = 0; i < n; ++i)
    {
        rhs[i] = b[i + 1] - b[0] * a[i + 1];
    }

    // 求解: zi = inv(sp) * rhs
    auto sp_inv = matrix_inverse(sp);
    return matrix_vector_multiply(sp_inv, rhs);
}

/**
 * @brief 使用初始条件应用滤波器（完全按照版本2的实现）
 */
inline std::vector<double> filter_with_zi(std::span<const double> signal,
                                          const FilterCoefficients &coeffs,
                                          std::vector<double> zi)
{
    if (coeffs.a.empty())
    {
        throw std::invalid_argument("Feedback filter coefficients are empty");
    }

    // 归一化系数
    std::vector<double> a = coeffs.a;
    std::vector<double> b = coeffs.b;

    double a0 = a[0];
    if (a0 == 0.0)
    {
        throw std::invalid_argument(
            "First feedback coefficient must be non-zero");
    }

    if (a0 != 1.0)
    {
        for (auto &val : a)
            val /= a0;
        for (auto &val : b)
            val /= a0;
    }

    const size_t input_size = signal.size();
    size_t filter_order = std::max(a.size(), b.size());

    a.resize(filter_order, 0.0);
    b.resize(filter_order, 0.0);
    zi.resize(filter_order, 0.0);

    std::vector<double> output(input_size);

    // 版本2的滤波实现
    for (size_t i = 0; i < input_size; ++i)
    {
        size_t order = filter_order - 1;
        while (order > 0)
        {
            if (i >= order)
            {
                zi[order - 1] = b[order] * signal[i - order]
                                - a[order] * output[i - order] + zi[order];
            }
            --order;
        }
        output[i] = b[0] * signal[i] + zi[0];
    }

    return output;
}

/**
 * @brief 零相位滤波（完全按照版本2的实现）
 */
inline std::vector<double> filtfilt(std::span<const double> signal,
                                    const FilterCoefficients &coeffs)
{
    const int len = static_cast<int>(signal.size());
    const int nfilt =
        static_cast<int>(std::max(coeffs.b.size(), coeffs.a.size()));
    const int nfact = 3 * (nfilt - 1);

    if (len <= nfact)
    {
        throw std::invalid_argument(
            "Input data too short! Must have length > 3 * filter_order");
    }

    // 计算初始条件
    auto zi_base = compute_filtfilt_zi(coeffs);

    // 左填充：2*signal[0] - signal[nfact:1:-1]
    std::vector<double> leftpad;
    for (int i = nfact; i >= 1; --i)
    {
        leftpad.push_back(2.0 * signal[0] - signal[i]);
    }

    // 右填充：2*signal[end] - signal[end-2:end-nfact-1:-1]
    std::vector<double> rightpad;
    for (int i = len - 2; i >= len - nfact - 1; --i)
    {
        rightpad.push_back(2.0 * signal[len - 1] - signal[i]);
    }

    // 组合信号
    std::vector<double> signal1;
    signal1.reserve(leftpad.size() + len + rightpad.size());
    signal1.insert(signal1.end(), leftpad.begin(), leftpad.end());
    signal1.insert(signal1.end(), signal.begin(), signal.end());
    signal1.insert(signal1.end(), rightpad.begin(), rightpad.end());

    // 前向滤波
    std::vector<double> zi = zi_base;
    double y0 = signal1[0];
    for (auto &z : zi)
    {
        z *= y0;
    }
    auto signal2 = filter_with_zi(signal1, coeffs, zi);

    // 反转
    std::reverse(signal2.begin(), signal2.end());

    // 后向滤波
    zi = zi_base;
    y0 = signal2[0];
    for (auto &z : zi)
    {
        z *= y0;
    }
    signal1 = filter_with_zi(signal2, coeffs, zi);

    // 提取结果（反转后的中间部分）
    std::vector<double> result;
    result.reserve(len);
    for (int i = signal1.size() - nfact - 1; i >= nfact; --i)
    {
        result.push_back(signal1[i]);
    }

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