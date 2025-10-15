#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>


#include "fft/fft.hpp"
#include "signal/butterworth_filter.hpp"
#include "signal/filter_apply.hpp"
#include "signal/filter_design.hpp"

using namespace msl;

std::vector<double> load_vector(const std::string &filename)
{
    std::ifstream fin(filename);
    std::vector<double> data;
    double val;
    while (fin >> val)
        data.push_back(val);
    return data;
}

std::vector<std::complex<double>> load_complex(const std::string &filename)
{
    std::ifstream fin(filename);
    std::vector<std::complex<double>> data;
    double re, im;
    while (fin >> re >> im)
        data.emplace_back(re, im);
    return data;
}

void save_vector(const std::string &filename, const std::vector<double> &data)
{
    std::ofstream fout(filename);
    for (auto v : data)
        fout << v << "\n";
}

void save_complex(const std::string &filename,
                  const std::vector<std::complex<double>> &data)
{
    std::ofstream fout(filename);
    for (auto c : data)
        fout << c.real() << " " << c.imag() << "\n";
}

int test_signal()
{
    // 读取信号
    std::ifstream fin("test_signal.txt");
    std::vector<double> t, x;
    double tval, xval;
    while (fin >> tval >> xval)
    {
        t.push_back(tval);
        x.push_back(xval);
    }

    // ========== FFT / IFFT ==========
    auto X_cpp = fft::fft(x);
    auto x_ifft_cpp = fft::ifft_real(X_cpp, x.size());
    save_complex("cpp_fft.txt", X_cpp);
    save_vector("cpp_ifft.txt", x_ifft_cpp);

    // ========== 滤波器设计 ==========
    // Test various configurations
    struct TestConfig
    {
        std::string name;
        int order;
        double fc_norm;
        signal::FilterType type;
    };

    /*std::vector<TestConfig> configs = {
        {"2nd order lowpass (fc=0.2)", 2, 0.2, signal::FilterType::lowpass},
        {"4th order lowpass (fc=0.2)", 4, 0.2, signal::FilterType::lowpass},
        {"4th order highpass (fc=0.3)", 4, 0.3, signal::FilterType::highpass},
    };
    for (const auto &config : configs)
    {
        signal::ButterworthFilter filt(
            config.order, config.fc_norm, config.type);
        auto coeffs = filt.coefficients();
        std::string fname = "cpp_" + config.name;
        std::replace(fname.begin(), fname.end(), ' ', '_');
        save_vector(fname + "_b.txt", coeffs.b);
        save_vector(fname + "_a.txt", coeffs.a);
    }*/

    // ========== 滤波器 ==========
    // 从 MATLAB 保存的系数读取
    /*std::ifstream fcoeff("filter_coeffs.txt");
    std::vector<double> b, a;
    double val;
    for (int i = 0; i < 5; ++i)
    {
        fcoeff >> val;
        b.push_back(val);
        fcoeff >> val;
        a.push_back(val);
    }
    msl::signal::FilterCoefficients coeffs{b, a};*/

    signal::ButterworthFilter filt(4, 0.05, 0.5, signal::FilterType::bandpass);
    // signal::ButterworthFilter filt(4, 0.02, signal::FilterType::highpass);
    auto coeffs = filt.coefficients();

    std::cout << "Filter order: " << coeffs.order() << std::endl;
    std::cout << "Numerator coefficients (b): ";
    for (auto v : coeffs.b)
        std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "Denominator coefficients (a): ";
    for (auto v : coeffs.a)
        std::cout << v << " ";
    std::cout << std::endl;

    auto y_filter_cpp = signal::filter(x, coeffs);
    auto y_filtfilt_cpp = signal::filtfilt(x, coeffs);
    save_vector("cpp_filter.txt", y_filter_cpp);
    save_vector("cpp_filtfilt.txt", y_filtfilt_cpp);

    std::cout << "Completed" << std::endl;
    return 0;
}
