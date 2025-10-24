#include <vector>

#include "matrix.hpp"
#include "signal.hpp"
#include "utils/data_io.hpp"

using namespace msl;
int test_signal()
{
    auto ori_data = utils::ReadData("KunmingSSJY.txt", 6, 3e4);
    auto data_1d =
        std::vector<double>(ori_data.begin(), ori_data.begin() + 3e4);
    matrixd ori_matrix(3e4, 6, std::span<const double>(ori_data));
    try
    {
        // test fourier domain bandpass filter
        auto fft_filter = signal::FourierDomainFilter(0.01, 0.8, 0);
        auto fft_filt = fft_filter.apply(data_1d);
        utils::WriteData(
            "test_result/signal/signal_fft_bandpass.txt", fft_filt, 1, 3e4);

        // test butterworth bandpass filter
        auto butter_filt = signal::butterworth_bandpass(data_1d, 4, 0.01, 0.8);
        utils::WriteData("test_result/signal/signal_butterworth_bandpass.txt",
                         butter_filt,
                         1,
                         3e4);

        // test matrix input with butterworth filter
        auto butter = signal::ButterworthFilter(4, 0.01, 0.8);
        matrixd butter_matrix_filt =
            signal::filtfilt_columns(ori_matrix, butter.coefficients());
        utils::WriteData(
            "test_result/signal/signal_butterworth_bandpass_matrix.txt",
            std::vector<double>(butter_matrix_filt.data(),
                                butter_matrix_filt.data()
                                    + butter_matrix_filt.size()),
            butter_matrix_filt.cols(),
            butter_matrix_filt.rows());
    }
    catch (const std::exception &e)
    {
        std::cerr << "Signal test failed: " << e.what() << std::endl;
        return -1;
    }

    std::cout << "Signal test completed successfully." << std::endl;
    return 0;
}


int main() { return test_signal(); }