/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: power_spectral_density.hpp
** -----
** File Created: Friday, 9th January 2026 14:58:10
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Friday, 6th March 2026 09:45:11
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_POWER_SPECTRAL_DENSITY_HPP
#define MSL_POWER_SPECTRAL_DENSITY_HPP

#include <algorithm>
#include <complex>
#include <span>
#include <stdexcept>
#include <vector>

#include "fft.hpp"
#include "window.hpp"

namespace msl::signal
{
// ========================================================================
// 1. Cross Power Spectral Density (CPSD) computation
// ========================================================================
/**
 * @brief Compute cross power spectral density using Welch's method.
 *
 * @param x First input signal.
 * @param y Second input signal.
 * @param output Output CPSD spectrum with size `nperseg`.
 * @param window Segment window coefficients, size must equal `nperseg`.
 * @param noverlap Overlap samples between adjacent segments.
 * @param nperseg Segment length and FFT length.
 */
inline void cpsd_welch(std::span<const double> x,
                       std::span<const double> y,
                       std::span<std::complex<double>> output,
                       const std::vector<double> &window = hann_window(1024),
                       size_t noverlap = 512,
                       size_t nperseg = 1024)
{
    if (nperseg == 0)
    {
        throw std::invalid_argument(
            "CPSD Welch: nperseg must be greater than 0");
    }
    if (noverlap >= nperseg)
    {
        throw std::invalid_argument(
            "CPSD Welch: noverlap must be less than nperseg");
    }
    if (x.size() != y.size())
    {
        throw std::invalid_argument("Input signals must have the same length.");
    }
    if (window.size() != nperseg)
    {
        throw std::invalid_argument(
            "CPSD Welch: window size must be equal to nperseg");
    }
    if (output.size() != nperseg)
    {
        throw std::invalid_argument(
            "CPSD Welch: output buffer size must be equal to nperseg");
    }
    if (x.size() < nperseg || y.size() < nperseg)
    {
        throw std::invalid_argument(
            "CPSD Welch: input signals must be at least as long "
            "as nperseg.");
    }
    size_t step = nperseg - noverlap;
    size_t num_segments = (std::min(x.size(), y.size()) - noverlap) / step;

    std::vector<std::complex<double>> psd_accum(nperseg,
                                                std::complex<double>(0.0, 0.0));
    double window_norm = 0.0;
    for (double w : window)
    {
        window_norm += w * w;
    }
    if (window_norm == 0.0)
    {
        throw std::invalid_argument(
            "CPSD Welch: window energy must be greater than 0");
    }

    for (size_t seg = 0; seg < num_segments; ++seg)
    {
        size_t start = seg * step;

        std::vector<double> x_segment(x.begin() + start,
                                      x.begin() + start + nperseg);
        std::vector<double> y_segment(y.begin() + start,
                                      y.begin() + start + nperseg);

        // Apply window
        for (size_t i = 0; i < nperseg; ++i)
        {
            x_segment[i] *= window[i];
            y_segment[i] *= window[i];
        }

        // Compute FFTs
        auto Xf = signal::fft(x_segment, nperseg);
        auto Yf = signal::fft(y_segment, nperseg);

        // Accumulate cross power
        for (size_t k = 0; k < nperseg; ++k)
        {
            psd_accum[k] += Xf[k] * std::conj(Yf[k]);
        }
    }

    // Average and normalize
    for (size_t k = 0; k < nperseg; ++k)
    {
        output[k] =
            psd_accum[k] / static_cast<double>(num_segments * window_norm);
    }
}

/**
 * @brief Return cross power spectral density using Welch's method.
 *
 * @param x First input signal.
 * @param y Second input signal.
 * @param window Segment window coefficients, size must equal `nperseg`.
 * @param noverlap Overlap samples between adjacent segments.
 * @param nperseg Segment length and FFT length.
 * @return CPSD spectrum with size `nperseg`.
 */
inline std::vector<std::complex<double>>
cpsd_welch(std::span<const double> x,
           std::span<const double> y,
           const std::vector<double> &window = hann_window(1024),
           size_t noverlap = 512,
           size_t nperseg = 1024)
{
    std::vector<std::complex<double>> output(nperseg);
    cpsd_welch(x, y, output, window, noverlap, nperseg);
    return output;
}

/**
 * @brief Compute cross power spectral density from one FFT frame.
 *
 * @param x First input signal.
 * @param y Second input signal.
 * @param output Output CPSD spectrum with size `nfft`.
 * @param nfft FFT length.
 */
inline void cpsd(std::span<const double> x,
                 std::span<const double> y,
                 std::span<std::complex<double>> output,
                 size_t nfft)
{
    if (nfft == 0)
    {
        throw std::invalid_argument("CPSD: nfft must be greater than 0");
    }
    if (x.size() != y.size())
    {
        throw std::invalid_argument("Input signals must have the same length.");
    }
    if (output.size() != nfft)
    {
        throw std::invalid_argument(
            "CPSD: output buffer size must be equal to nfft");
    }
    if (x.size() < nfft || y.size() < nfft)
    {
        throw std::invalid_argument(
            "CPSD: input signals must be at least as long as nfft.");
    }

    // Compute FFTs
    auto Xf = signal::fft(x, nfft);
    auto Yf = signal::fft(y, nfft);

    // Compute Cross Power Spectral Density
    for (size_t k = 0; k < nfft; ++k)
    {
        output[k] = Xf[k] * std::conj(Yf[k]);
    }
}

/**
 * @brief Return cross power spectral density from one FFT frame.
 *
 * @param x First input signal.
 * @param y Second input signal.
 * @param nfft FFT length.
 * @return CPSD spectrum with size `nfft`.
 */
inline std::vector<std::complex<double>>
cpsd(std::span<const double> x, std::span<const double> y, size_t nfft)
{
    std::vector<std::complex<double>> output(nfft);
    cpsd(x, y, output, nfft);
    return output;
}

// ========================================================================
// 2. Power Spectral Density (PSD) computation
// ========================================================================
/**
 * @brief Compute power spectral density using Welch's method.
 *
 * @param x Input signal.
 * @param output Output PSD spectrum with size `nperseg`.
 * @param window Segment window coefficients, size must equal `nperseg`.
 * @param noverlap Overlap samples between adjacent segments.
 * @param nperseg Segment length and FFT length.
 */
inline void psd_welch(std::span<const double> x,
                      std::span<double> output,
                      const std::vector<double> &window = hann_window(1024),
                      size_t noverlap = 512,
                      size_t nperseg = 1024)
{
    if (output.size() != nperseg)
    {
        throw std::invalid_argument(
            "PSD Welch: output buffer size must be equal to nperseg");
    }

    std::vector<std::complex<double>> cpsd_result(nperseg);
    cpsd_welch(x, x, cpsd_result, window, noverlap, nperseg);

    for (size_t k = 0; k < nperseg; ++k)
    {
        output[k] = std::real(cpsd_result[k]);
    }
}

/**
 * @brief Return power spectral density using Welch's method.
 *
 * @param x Input signal.
 * @param window Segment window coefficients, size must equal `nperseg`.
 * @param noverlap Overlap samples between adjacent segments.
 * @param nperseg Segment length and FFT length.
 * @return PSD spectrum with size `nperseg`.
 */
inline std::vector<double>
psd_welch(std::span<const double> x,
          const std::vector<double> &window = hann_window(1024),
          size_t noverlap = 512,
          size_t nperseg = 1024)
{
    std::vector<double> output(nperseg);
    psd_welch(x, output, window, noverlap, nperseg);
    return output;
}

/**
 * @brief Compute power spectral density from one FFT frame.
 *
 * @param x Input signal.
 * @param output Output PSD spectrum with size `nfft`.
 * @param nfft FFT length.
 */
inline void
psd(std::span<const double> x, std::span<double> output, size_t nfft)
{
    if (nfft == 0)
    {
        throw std::invalid_argument("PSD: nfft must be greater than 0");
    }
    if (output.size() != nfft)
    {
        throw std::invalid_argument(
            "PSD: output buffer size must be equal to nfft");
    }

    std::vector<std::complex<double>> cpsd_result(nfft);
    cpsd(x, x, cpsd_result, nfft);

    for (size_t k = 0; k < nfft; ++k)
    {
        output[k] = std::real(cpsd_result[k]);
    }
}

/**
 * @brief Return power spectral density from one FFT frame.
 *
 * @param x Input signal.
 * @param nfft FFT length.
 * @return PSD spectrum with size `nfft`.
 */
inline std::vector<double> psd(std::span<const double> x, size_t nfft)
{
    std::vector<double> output(nfft);
    psd(x, output, nfft);
    return output;
}

} // namespace msl::signal

#endif // MSL_POWER_SPECTRAL_DENSITY_HPP