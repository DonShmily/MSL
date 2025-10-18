#ifndef MSL_TEST_HELPER_HPP
#define MSL_TEST_HELPER_HPP

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>


#include "fft/fft.hpp" // For frequency analysis
#include "signal/filter_apply.hpp"
#include "signal/filter_design.hpp"

using namespace msl;
using namespace msl::signal;

void print_coefficients(const FilterCoefficients &coeffs,
                        const std::string &name = "Filter");
void print_vector(const std::vector<double> &v,
                  const std::string &name,
                  size_t max_print = 10);
void print_complex_vector(const std::vector<std::complex<double>> &v,
                          const std::string &name = "Complex Vector",
                          size_t max_print = 10);
std::vector<double> read_csv(const std::string &filename);
double max_abs_error(const std::vector<double> &a,
                     const std::vector<double> &b);
double max_abs_error_complex(const std::vector<std::complex<double>> &a,
                             const std::vector<std::complex<double>> &b);
double rms_error(const std::vector<double> &a, const std::vector<double> &b);
void print_result(const std::string &test_name, double error, double threshold);

// 从文件中读取矩阵数据
// @param file_name 文件名
// @return 读取到的矩阵数据
std::vector<double> ReadData(const std::string &file_name, int col, int row);

// 将矩阵数据写入文件
// @param file_name 文件名
// @param data 矩阵数据
void WriteData(const std::string &file_name,
               const std::vector<double> &data,
               int col,
               int row);

// 读取JSON至字符串
// @param file_name 文件名
std::string ReadJson(const std::string &file_name);

#endif // MSL_TEST_HELPER_HPP