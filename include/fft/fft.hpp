/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025
**  File: fft.hpp
**  Description: FFT interface using Eigen's unsupported FFT module
*/

#ifndef MSL_FFT_HPP
#define MSL_FFT_HPP

#define _USE_MATH_DEFINES

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <cassert>
#include <complex>
#include <cstddef>

#include <span>
#include <vector>

#include <eigen3/unsupported/Eigen/FFT>

#include "matrix/complex_matrix_owned.hpp"
#include "matrix/real_matrix_owned.hpp"

namespace msl::fft
{

// ============================================================================
// 1. 1D FFT for vectors
// ============================================================================

/**
 * @brief 1D Forward FFT: real vector -> complex vector
 */
inline std::vector<std::complex<double>> fft(std::span<const double> input,
                                             std::size_t nfft = 0)
{
    nfft = (nfft == 0) ? input.size() : nfft;
    Eigen::FFT<double> fft_engine;
    std::vector<std::complex<double>> output(nfft);

    // Perform FFT
    fft_engine.fwd(output.data(), input.data(), nfft);

    return output;
}

/**
 * @brief 1D Forward FFT: real vector -> complex vector (vector overload)
 */
inline std::vector<std::complex<double>> fft(const std::vector<double> &input,
                                             std::size_t nfft = 0)
{
    return fft(std::span<const double>(input), nfft);
}

/**
 * @brief 1D Forward FFT: complex vector -> complex vector
 *
 * Full complex-to-complex FFT (no symmetry assumptions)
 *
 * @param input Complex-valued input signal
 * @return Complex-valued frequency domain (same size as input)
 */
inline std::vector<std::complex<double>>
fft(std::span<const std::complex<double>> input, std::size_t nfft = 0)
{
    nfft = (nfft == 0) ? input.size() : nfft;

    Eigen::FFT<double> fft_engine;
    std::vector<std::complex<double>> output(nfft);

    fft_engine.fwd(output.data(), input.data(), nfft);

    return output;
}

/**
 * @brief 1D Forward FFT: complex vector -> complex vector (vector overload)
 */
inline std::vector<std::complex<double>>
fft(const std::vector<std::complex<double>> &input, std::size_t nfft = 0)
{
    return fft(std::span<const std::complex<double>>(input), nfft);
}

// ============================================================================
// 2. 1D Inverse FFT for vectors
// ============================================================================

/**
 * @brief 1D Inverse FFT: complex vector -> complex vector
 *
 * @param input Complex frequency domain data
 * @return Complex time domain signal
 */
inline std::vector<std::complex<double>>
ifft(std::span<const std::complex<double>> input, std::size_t nfft = 0)
{
    nfft = (nfft == 0) ? input.size() : nfft;
    Eigen::FFT<double> fft_engine;
    std::vector<std::complex<double>> output(nfft);

    fft_engine.inv(output.data(), input.data(), nfft);

    return output;
}

/**
 * @brief 1D Inverse FFT: complex vector -> complex vector (vector overload)
 */
inline std::vector<std::complex<double>>
ifft(const std::vector<std::complex<double>> &input, std::size_t nfft = 0)
{
    return ifft(std::span<const std::complex<double>>(input), nfft);
}

/**
 * @brief 1D Inverse FFT: complex vector -> real vector
 *
 * Assumes input has conjugate symmetry (from real FFT)
 *
 *
 * @warning Input must have proper conjugate symmetry for real output
 */
inline std::vector<double>
ifft_real(std::span<const std::complex<double>> input, size_t nfft)
{
    Eigen::FFT<double> fft_engine;
    std::vector<double> output(nfft);

    fft_engine.inv(output.data(), input.data(), nfft);

    return output;
}

/**
 * @brief 1D Inverse FFT: complex vector -> real vector (vector overload)
 */
inline std::vector<double>
ifft_real(const std::vector<std::complex<double>> &input, size_t nfft)
{
    return ifft_real(std::span<const std::complex<double>>(input), nfft);
}

// ============================================================================
// 3. 2D FFT for matrices (column-wise FFT)
// ============================================================================

/**
 * @brief 2D FFT: real matrix -> complex matrix (column-wise)
 *
 * Applies 1D FFT to each column of the matrix independently.
 * Useful for processing multiple signals simultaneously.
 */
inline matrixc fft_columns(const matrixd &input, size_t nfft = 0)
{
    size_t n_rows = input.rows();
    size_t n_cols = input.cols();

    nfft = (nfft == 0) ? n_rows : nfft;

    matrixc output(nfft, n_cols);

    Eigen::FFT<double> fft_engine;

    // Process each column
    for (size_t j = 0; j < n_cols; ++j)
    {
        auto col_span = input.column(j);

        auto col_fft = output.column(j);
        fft_engine.fwd(col_fft.data(), col_span.data(), nfft);
    }

    return output;
}

/**
 * @brief 2D FFT: complex matrix -> complex matrix (column-wise)
 *
 * Full complex FFT for each column
 */
inline matrixc fft_columns(const matrixc &input, size_t nfft = 0)
{
    size_t n_rows = input.rows();
    size_t n_cols = input.cols();

    nfft = (nfft == 0) ? n_rows : nfft;

    matrixc output(nfft, n_cols);

    Eigen::FFT<double> fft_engine;

    for (size_t j = 0; j < n_cols; ++j)
    {
        auto col_span = input.column(j);
        auto col_fft = output.column(j);
        fft_engine.fwd(col_fft.data(), col_span.data(), nfft);
    }

    return output;
}

// ============================================================================
// 4. 2D Inverse FFT for matrices (column-wise)
// ============================================================================

/**
 * @brief 2D Inverse FFT: complex matrix -> complex matrix (column-wise)
 */
inline matrixc ifft_columns(const matrixc &input, size_t nfft = 0)
{
    size_t n_rows = input.rows();
    size_t n_cols = input.cols();

    nfft = (nfft == 0) ? n_rows : nfft;

    matrixc output(nfft, n_cols);

    Eigen::FFT<double> fft_engine;

    for (size_t j = 0; j < n_cols; ++j)
    {
        auto col_span = input.column(j);
        auto col_ifft = output.column(j);
        fft_engine.inv(col_ifft.data(), col_span.data(), nfft);
    }

    return output;
}

/**
 * @brief 2D Inverse FFT: complex matrix -> real matrix (column-wise)
 *
 * @param input Complex frequency domain matrix
 * @param output_rows Expected number of rows in output (original signal length)
 */
inline matrixd ifft_columns_real(const matrixc &input, size_t nfft)
{
    size_t n_cols = input.cols();

    nfft = (nfft == 0) ? n_cols : nfft;

    matrixd output(nfft, n_cols);

    Eigen::FFT<double> fft_engine;

    for (size_t j = 0; j < n_cols; ++j)
    {
        auto col_span = input.column(j);
        auto col_ifft = output.column(j);
        fft_engine.inv(col_ifft.data(), col_span.data(), nfft);
    }

    return output;
}

// ============================================================================
// 5. 2D FFT for matrices (row-wise FFT)
// ============================================================================

/**
 * @brief 2D FFT: real matrix -> complex matrix (row-wise)
 *
 * Applies 1D FFT to each row of the matrix independently.
 *
 * @param input Real-valued matrix (each row is a signal)
 * @return Complex matrix with FFT of each row
 */
inline matrixc fft_rows(const matrixd &input, size_t nfft = 0)
{
    size_t n_rows = input.rows();
    size_t n_cols = input.cols();

    nfft = (nfft == 0) ? n_cols : nfft;

    matrixc output(n_rows, nfft);

    Eigen::FFT<double> fft_engine;

    for (size_t i = 0; i < n_rows; ++i)
    {
        auto row_vec = input.get_row(i);
        auto row_fft = output.get_row(i);
        fft_engine.fwd(row_fft.data(), row_vec.data(), nfft);
    }

    return output;
}

/**
 * @brief 2D Inverse FFT: complex matrix -> real matrix (row-wise)
 */
inline matrixd ifft_rows_real(const matrixc &input, size_t nfft)
{
    size_t n_rows = input.rows();

    matrixd output(n_rows, nfft);

    Eigen::FFT<double> fft_engine;

    for (size_t i = 0; i < n_rows; ++i)
    {
        auto row_vec = input.get_row(i);
        auto row_ifft = output.get_row(i);
        fft_engine.inv(row_ifft.data(), row_vec.data(), nfft);
    }

    return output;
}

// ============================================================================
// 6. Utility functions
// ============================================================================

/**
 * @brief Compute FFT frequencies for a given sample rate
 *
 * @param n Number of samples
 * @param sample_rate Sampling rate (Hz)
 * @return Vector of frequency values
 */
inline std::vector<double> fft_frequencies(size_t n, double sample_rate = 1.0)
{
    std::vector<double> freqs(n);
    double df = sample_rate / n;

    for (size_t i = 0; i < freqs.size(); ++i)
    {
        freqs[i] = i * df;
    }

    return freqs;
}

/**
 * @brief Compute power spectrum from FFT result
 *
 * @param fft_result Complex FFT output
 * @return Power spectrum (|FFT|^2)
 */
inline std::vector<double>
power_spectrum(std::span<const std::complex<double>> fft_result)
{
    std::vector<double> power(fft_result.size());

    for (size_t i = 0; i < fft_result.size(); ++i)
    {
        power[i] = std::norm(fft_result[i]); // |z|^2
    }

    return power;
}

/**
 * @brief Compute magnitude spectrum from FFT result
 *
 * @param fft_result Complex FFT output
 * @return Magnitude spectrum (|FFT|)
 */
inline std::vector<double>
magnitude_spectrum(std::span<const std::complex<double>> fft_result)
{
    std::vector<double> magnitude(fft_result.size());

    for (size_t i = 0; i < fft_result.size(); ++i)
    {
        magnitude[i] = std::abs(fft_result[i]);
    }

    return magnitude;
}

/**
 * @brief Compute phase spectrum from FFT result
 *
 * @param fft_result Complex FFT output
 * @return Phase spectrum (arg(FFT))
 */
inline std::vector<double>
phase_spectrum(std::span<const std::complex<double>> fft_result)
{
    std::vector<double> phase(fft_result.size());

    for (size_t i = 0; i < fft_result.size(); ++i)
    {
        phase[i] = std::arg(fft_result[i]);
    }

    return phase;
}

/**
 * @brief Apply window function to signal (Hamming window)
 *
 * @param signal Input signal (will be modified in-place)
 */
inline void apply_hamming_window(std::span<double> signal)
{
    size_t n = signal.size();

    for (size_t i = 0; i < n; ++i)
    {
        double window = 0.54 - 0.46 * std::cos(2.0 * M_PI * i / (n - 1));
        signal[i] *= window;
    }
}

/**
 * @brief Apply window function to signal (Hann window)
 *
 * @param signal Input signal (will be modified in-place)
 */
inline void apply_hann_window(std::span<double> signal)
{
    size_t n = signal.size();

    for (size_t i = 0; i < n; ++i)
    {
        double window = 0.5 * (1.0 - std::cos(2.0 * M_PI * i / (n - 1)));
        signal[i] *= window;
    }
}

} // namespace msl::fft

#endif // MSL_FFT_HPP