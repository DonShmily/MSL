#include <complex>
#include <cstddef>
#include <vector>
#include "helper.hpp"
#include "test.hpp"

#include "fft.hpp"
#include "signal/fourier_domain_filter.hpp"

void test_fft_filter()
{
    int row = 30000;
    int col = 1;
    auto data = ReadData("sig.txt", col, row);

    msl::signal::FourierDomainFilter filter(
        0.1, 0.3, 0, msl::signal::FilterType::bandpass);
    filter.set_window_type(
        msl::signal::FourierDomainFilter::WindowType::rectangular);
    filter.set_transition_band(0.000);
    auto filtered_data = filter.apply(data);
    WriteData("fft_filter.txt", filtered_data, col, row);

    auto fft_data = fft::fft(data);
    std::vector<std::complex<double>> filtered_fft2(fft_data.size());
    std::size_t low_idx = static_cast<std::size_t>(0.1 * row / 2);
    std::size_t high_idx = static_cast<std::size_t>(0.3 * row / 2);
    for (std::size_t i = 0; i < filtered_fft2.size(); ++i)
    {
        if (i < low_idx || i > high_idx && i < row - high_idx
            || i > row - low_idx)
        {
            filtered_fft2[i] = 0;
        }
        else
        {
            filtered_fft2[i] = fft_data[i];
        }
    }
    auto ifft_data2 = fft::ifft_real(filtered_fft2);
    WriteData("fft_filter_manual_freq_response.txt", ifft_data2, col, row);
}
