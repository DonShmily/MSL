/**
**            qREST - Quick Response Evaluation for Safety Tagging
**     Institute of Engineering Mechanics, China Earthquake Administration
**
**                 Copyright 2024 - 2024 QLab, Dong Feiyue
**                          All Rights Reserved.
**
** Project: qREST
** File: \src\numerical_algorithm\butterworth_filter_design.cpp
** -----
** File Created: Friday, 12th July 2024 22:17:18
** Author: Dong Feiyue (donfeiyue@outlook.com)
** -----
** Last Modified: Saturday, 10th August 2024 22:28:11
** Modified By: Dong Feiyue (donfeiyue@outlook.com)
*/

// Description:
// 数值算法：Butterworth滤波器设计类的实现


// associated header
#include "butterworth_filter_design.h"

// stdc++ headers
#include <cmath>
#include <complex>
#include <vector>

// project headers
#include "basic_filter_design.h"


namespace numerical_algorithm
{

// 通过滤波器阶数和滤波器类型构造
ButterworthFilterDesign::ButterworthFilterDesign(int filter_order,
                                                 FilterType filter_type)
    : BasicFilterDesign(filter_order, 0, 0, filter_type)
{
    DesignFilter();
}

// 带通滤波器构造函数
ButterworthFilterDesign::ButterworthFilterDesign(int filter_order,
                                                 double low_frequency,
                                                 double high_frequency,
                                                 FilterType filter_type)
    : BasicFilterDesign(filter_order,
                        low_frequency,
                        high_frequency,
                        filter_type)
{
    DesignFilter();
}

// 低通/高通滤波器构造函数（默认低通）
ButterworthFilterDesign::ButterworthFilterDesign(int filter_order,
                                                 double frequency,
                                                 FilterType filter_type)
    : BasicFilterDesign(filter_order, frequency, filter_type)
{
    DesignFilter();
}

// 计算滤波器参数的函数
void ButterworthFilterDesign::ComputeDenCoeffs()
{
    double cp = std::cos(M_PI * (high_frequency_ + low_frequency_)
                         / 2.0); // cosine of phi
    double theta = M_PI * (high_frequency_ - low_frequency_)
                   / 2.0; // M_PI * (high_frequency_ - low_frequency_) / 2.0
    double st = std::sin(theta);                    // sine of theta
    double ct = std::cos(theta);                    // cosine of theta
    double s2t = 2.0 * st * ct;                     // sine of 2*theta
    double c2t = 2.0 * ct * ct - 1.0;               // cosine of 2*theta
    std::vector<double> rCoeffs(2 * filter_order_); // z^-2 coefficients
    std::vector<double> tCoeffs(2 * filter_order_); // z^-1 coefficients
    double poleAngle;                               // pole angle
    double sinPoleAngle;                            // sine of pole angle
    double cosPoleAngle;                            // cosine of pole angle
    double a;                                       // workspace variables

    for (int k = 0; k < filter_order_; ++k)
    {
        poleAngle = M_PI * (double)(2 * k + 1) / (double)(2 * filter_order_);
        sinPoleAngle = std::sin(poleAngle);
        cosPoleAngle = std::cos(poleAngle);
        a = 1.0 + s2t * sinPoleAngle;
        rCoeffs[2 * k] = c2t / a;
        rCoeffs[2 * k + 1] = s2t * cosPoleAngle / a;
        tCoeffs[2 * k] = -2.0 * cp * (ct + st * sinPoleAngle) / a;
        tCoeffs[2 * k + 1] = -2.0 * cp * st * cosPoleAngle / a;
    }

    coefficients_a_ = TrinomialMultiply(tCoeffs, rCoeffs);

    coefficients_a_[1] = coefficients_a_[0];
    coefficients_a_[0] = 1.0;
    for (std::size_t k = 3; k <= 2 * filter_order_; ++k)
        coefficients_a_[k] = coefficients_a_[2 * k - 2];

    for (auto i = coefficients_a_.size(); i > filter_order_ * 2 + 1; i--)
        coefficients_a_.pop_back();
}

void ButterworthFilterDesign::ComputeNumCoeffs()
{
    std::vector<double> tCoeffs;
    coefficients_b_.resize(2 * filter_order_ + 1, 0);
    std::vector<std::complex<double>> normalizedKernel(2 * filter_order_ + 1);

    std::vector<double> numbers;
    for (double n = 0; n < filter_order_ * 2 + 1; n++) numbers.push_back(n);
    int i;

    tCoeffs = ComputeHP();

    for (i = 0; i < filter_order_; ++i)
    {
        coefficients_b_[2 * i] = tCoeffs[i];
        coefficients_b_[2 * i + 1] = 0.0;
    }
    coefficients_b_[2 * filter_order_] = tCoeffs[filter_order_];

    double cp[2] = {};
    double wn;
    cp[0] = 2 * 2.0 * tan(M_PI * low_frequency_ / 2.0);
    cp[1] = 2 * 2.0 * tan(M_PI * high_frequency_ / 2.0);

    // center frequency
    wn = std::sqrt(cp[0] * cp[1]);
    wn = 2 * std::atan2(wn, 4);
    const std::complex<double> result = std::complex<double>(-1, 0);

    for (int k = 0; k < filter_order_ * 2 + 1; k++)
    {
        normalizedKernel[k] = std::exp(-sqrt(result) * wn * numbers[k]);
    }
    double b = 0;
    double den = 0;
    for (int d = 0; d < filter_order_ * 2 + 1; d++)
    {
        b += std::real(normalizedKernel[d] * coefficients_b_[d]);
        den += std::real(normalizedKernel[d] * coefficients_a_[d]);
    }
    for (int c = 0; c < filter_order_ * 2 + 1; c++)
    {
        coefficients_b_[c] = (coefficients_b_[c] * den) / b;
    }

    for (auto i = coefficients_b_.size() - 1; i > filter_order_ * 2 + 1; i--)
        coefficients_b_.pop_back();
}

std::vector<double>
ButterworthFilterDesign::TrinomialMultiply(const std::vector<double> &b,
                                           const std::vector<double> &c) const
{
    int i, j;
    std::vector<double> retVal(4 * filter_order_);

    retVal[0] = b[0];
    retVal[1] = b[1];
    retVal[2] = c[0];
    retVal[3] = c[1];

    for (i = 1; i < filter_order_; ++i)
    {
        retVal[2 * (2 * i + 1)] += c[2 * i] * retVal[2 * (2 * i - 1)]
                                   - c[2 * i + 1] * retVal[2 * (2 * i - 1) + 1];
        retVal[2 * (2 * i + 1) + 1] += c[2 * i] * retVal[2 * (2 * i - 1) + 1]
                                       + c[2 * i + 1] * retVal[2 * (2 * i - 1)];

        for (j = 2 * i; j > 1; --j)
        {
            retVal[2 * j] += b[2 * i] * retVal[2 * (j - 1)]
                             - b[2 * i + 1] * retVal[2 * (j - 1) + 1]
                             + c[2 * i] * retVal[2 * (j - 2)]
                             - c[2 * i + 1] * retVal[2 * (j - 2) + 1];
            retVal[2 * j + 1] += b[2 * i] * retVal[2 * (j - 1) + 1]
                                 + b[2 * i + 1] * retVal[2 * (j - 1)]
                                 + c[2 * i] * retVal[2 * (j - 2) + 1]
                                 + c[2 * i + 1] * retVal[2 * (j - 2)];
        }

        retVal[2] += b[2 * i] * retVal[0] - b[2 * i + 1] * retVal[1] + c[2 * i];
        retVal[3] +=
            b[2 * i] * retVal[1] + b[2 * i + 1] * retVal[0] + c[2 * i + 1];
        retVal[0] += b[2 * i];
        retVal[1] += b[2 * i + 1];
    }

    return retVal;
}

std::vector<double> ButterworthFilterDesign::ComputeLP() const
{
    std::vector<double> numCoeffs(filter_order_ + 1);
    int m;
    int i;

    numCoeffs[0] = 1;
    numCoeffs[1] = filter_order_;
    m = filter_order_ / 2;
    for (i = 2; i <= m; ++i)
    {
        numCoeffs[i] = (double)(filter_order_ - i + 1) * numCoeffs[i - 1] / i;
        numCoeffs[filter_order_ - i] = numCoeffs[i];
    }
    numCoeffs[filter_order_ - 1] = filter_order_;
    numCoeffs[filter_order_] = 1;

    return numCoeffs;
}

std::vector<double> ButterworthFilterDesign::ComputeHP() const
{
    std::vector<double> numCoeffs;
    int i;

    numCoeffs = ComputeLP();

    for (i = 0; i <= filter_order_; ++i)
        if (i % 2) numCoeffs[i] = -numCoeffs[i];

    return numCoeffs;
}

} // namespace numerical_algorithm