#include <vector>

#include "matrix.hpp"
#include "signal.hpp"
#include "signal/detrend.hpp"
#include "utils/data_io.hpp"

using namespace msl;

int test_filter()
{
    auto ori_data = utils::ReadData("KunmingSSJY.txt", 6, 3e4);
    auto data_1d =
        std::vector<double>(ori_data.begin(), ori_data.begin() + 3e4);
    matrix::matrixd ori_matrix(3e4, 6, std::span<const double>(ori_data));
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
        matrix::matrixd butter_matrix_filt =
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

int test_fft()
{
    auto ori_data = utils::ReadData("KunmingSSJY.txt", 6, 3e4);
    auto data_1d =
        std::vector<double>(ori_data.begin(), ori_data.begin() + 3e4);
    matrix::matrixd ori_matrix(3e4, 6, std::span<const double>(ori_data));

    try
    {

        // test r2c fft
        auto fft_result = signal::fft(data_1d);
        utils::WriteComplexData(
            "test_result/signal/fft_result.txt", fft_result, 1, 3e4);

        // test c2c fft
        auto fft_result_c2c = signal::fft(fft_result);
        utils::WriteComplexData(
            "test_result/signal/fft_result_c2c.txt", fft_result_c2c, 1, 3e4);

        // test c2r ifft
        auto ifft_result = signal::ifft_real(fft_result, 3e4);
        utils::WriteData(
            "test_result/signal/ifft_result.txt", ifft_result, 1, 3e4);

        // test c2c ifft
        auto ifft_result_c2c = signal::ifft(fft_result_c2c, 3e4);
        utils::WriteComplexData(
            "test_result/signal/ifft_result_c2c.txt", ifft_result_c2c, 1, 3e4);

        // test matrix fft
        auto fft_matrix_result = signal::fft_columns(ori_matrix);
        auto res_vec = std::vector<std::complex<double>>(
            fft_matrix_result.data(),
            fft_matrix_result.data() + fft_matrix_result.size());
        utils::WriteComplexData("test_result/signal/fft_matrix_result.txt",
                                res_vec,
                                fft_matrix_result.cols(),
                                fft_matrix_result.rows());
    }
    catch (const std::exception &e)
    {
        std::cerr << "FFT test failed: " << e.what() << std::endl;
        return -1;
    }

    std::cout << "FFT test completed successfully." << std::endl;
    return 0;
}

int test_psd()
{
    auto ori_data = utils::ReadData("KunmingSSJY.txt", 6, 3e4);
    auto data_1d =
        std::vector<double>(ori_data.begin(), ori_data.begin() + 3e4);
    auto data_2d =
        std::vector<double>(ori_data.begin() + 3e4, ori_data.begin() + 6e4);
    try
    {
        // test psd
        auto psd_1d = signal::psd_welch(data_1d);
        utils::WriteData(
            "test_result/signal/psd_1d.txt", psd_1d, 1, psd_1d.size());

        // test cpsd
        auto cpsd_2d = signal::cpsd_welch(data_1d, data_2d);
        auto cpsd_2d_mag = std::vector<double>(cpsd_2d.size());
        for (size_t i = 0; i < cpsd_2d.size(); i++)
        {
            cpsd_2d_mag[i] = std::abs(cpsd_2d[i]);
        }
        // utils::WriteComplexData(
        // "test_result/signal/cpsd_2d.txt", cpsd_2d, 1, cpsd_2d.size());
        utils::WriteData("test_result/signal/cpsd_2d_mag.txt",
                         cpsd_2d_mag,
                         1,
                         cpsd_2d_mag.size());
    }
    catch (const std::exception &e)
    {
        std::cerr << "PSD test failed: " << e.what() << std::endl;
        return -1;
    }

    std::cout << "PSD test completed successfully." << std::endl;
    return 0;
}

int test_detrend()
{

    auto ori_data = utils::ReadData("err_data.txt", 1, 3e4);
    auto detrended = signal::detrend(ori_data, 2);

    try
    {
        utils::WriteData("test_result/signal/detrended_data.txt",
                         detrended,
                         1,
                         detrended.size());
    }
    catch (const std::exception &e)
    {
        std::cerr << "Detrend test failed: " << e.what() << std::endl;
        return -1;
    }

    std::cout << "Detrend test completed successfully." << std::endl;
    return 0;
}

int main()
{
    // return test_filter() + test_fft() + test_psd();
    return test_detrend();
}
