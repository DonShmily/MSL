/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: test_numerical_old.cpp
** -----
** File Created: Thursday, 8th January 2026 20:46:14
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Thursday, 8th January 2026 20:46:17
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/


#include "../numerical_old/filtfilt.h"
#include "utils/data_io.hpp"

int test_filtfilt()
{
    double fs = 50;
    double low = 0.1 / (fs / 2);
    double high = 10.0 / (fs / 2);
    // 生成滤波器生成器
    auto filter_generator = numerical_algorithm::ButterworthFilterDesign(
        4, low, high, numerical_algorithm::FilterType::bandpass);

    // 生成filtfilt滤波器
    auto filtfilt_filter = numerical_algorithm::FiltFilt(filter_generator);

    auto ori_data = msl::utils::ReadData("KunmingSSJY.txt", 6, 3e4);
    auto data_1d =
        std::vector<double>(ori_data.begin(), ori_data.begin() + 3e4);

    auto filtered_data = filtfilt_filter.Filtering(data_1d);

    // 保存滤波结果
    msl::utils::WriteData(
        "test_result/numerical_old/filtfilt_filtered_data.txt",
        filtered_data,
        1,
        filtered_data.size());

    return 0;
}

int main()
{
    test_filtfilt();
    return 0;
}