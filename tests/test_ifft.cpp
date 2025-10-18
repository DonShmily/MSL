#include <complex>
#include <vector>
#include "helper.hpp"
#include "test.hpp"

#include "fft.hpp"

void fft_ifft()
{
    int row = 30000;
    int col = 1;
    auto data = ReadData("sig.txt", col, row);

    auto fft_data = fft::fft(data);

    auto ifft_data = fft::ifft_real(fft_data);
    WriteData("ifft_real.txt", ifft_data, col, row);

    auto ifft_data2 = fft::ifft(fft_data);
    std::vector<double> ifft_data3(ifft_data2.size());
    for (size_t i = 0; i < ifft_data2.size(); ++i)
    {
        ifft_data3[i] = ifft_data2[i].real();
    }
    WriteData("ifft.txt", ifft_data3, col, row);

    std::vector<std::complex<double>> fft_half(fft_data.begin(),
                                               fft_data.end());
    for (auto it = fft_half.begin() + row / 2 + 1; it != fft_half.end(); ++it)
    {
        *it = 0;
    }
    auto ifft_data4 = fft::ifft_real(fft_half);
    WriteData("ifft_real_half.txt", ifft_data4, col, row);

    auto ifft_data5 = fft::ifft(fft_half);
    std::vector<double> ifft_data6(ifft_data5.size());
    for (size_t i = 0; i < ifft_data5.size(); ++i)
    {
        ifft_data6[i] = ifft_data5[i].real();
    }
    WriteData("ifft_half.txt", ifft_data6, col, row);
}

void test_ifft() { fft_ifft(); }