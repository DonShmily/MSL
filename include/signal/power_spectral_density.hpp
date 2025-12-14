/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: power_spectral_density.hpp
** -----
** File Created: Friday, 7th November 2025 20:44:37
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:05:38
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_POWER_SPECTRAL_DENSITY_HPP
#define MSL_POWER_SPECTRAL_DENSITY_HPP

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
// Compute Cross Power Spectral Density using Welch's method
inline std::vector<std::complex<double>>
cpsd_welch(std::span<const double> x,
           std::span<const double> y,
           const std::vector<double> &window = hann_window(1024),
           size_t noverlap = 512,
           size_t nperseg = 1024)
{
    if (x.size() < nperseg || y.size() < nperseg)
    {
        throw std::invalid_argument("Input signals must be at least as long "
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
        psd_accum[k] /= static_cast<double>(num_segments * window_norm);
    }

    return psd_accum;
}

// Compute Cross Power Spectral Density without Welch's method
inline std::vector<std::complex<double>>
cpsd(std::span<const double> x, std::span<const double> y, size_t nfft)
{
    if (x.size() != y.size())
    {
        throw std::invalid_argument("Input signals must have the same length.");
    }
    std::vector<double> x_padded(nfft, 0.0);
    std::vector<double> y_padded(nfft, 0.0);

    if (nfft < x.size())
    {
        std::copy(x.begin(), x.begin() + nfft, x_padded.begin());
        std::copy(y.begin(), y.begin() + nfft, y_padded.begin());
    }
    else
    {
        std::copy(x.begin(), x.end(), x_padded.begin());
        std::copy(y.begin(), y.end(), y_padded.begin());
    }

    // Compute FFTs
    auto Xf = signal::fft(x_padded, nfft);
    auto Yf = signal::fft(y_padded, nfft);

    // Compute Cross Power Spectral Density
    std::vector<std::complex<double>> cpsd(nfft);
    for (size_t k = 0; k < nfft; ++k)
    {
        cpsd[k] = Xf[k] * std::conj(Yf[k]);
    }

    return cpsd;
}

// ========================================================================
// 2. Power Spectral Density (PSD) computation
// ========================================================================
// Compute Power Spectral Density
inline std::vector<double>
psd_welch(std::span<const double> x,
          const std::vector<double> &window = hann_window(1024),
          size_t noverlap = 512,
          size_t nperseg = 1024)
{
    auto cpsd_result = cpsd_welch(x, x, window, noverlap, nperseg);
    std::vector<double> psd_result(cpsd_result.size());

    for (size_t k = 0; k < cpsd_result.size(); ++k)
    {
        psd_result[k] = std::real(cpsd_result[k]);
    }

    return psd_result;
}

// Compute Power Spectral Density without Welch's method
inline std::vector<double> psd(std::span<const double> x, size_t nfft)
{
    auto cpsd_result = cpsd(x, x, nfft);
    std::vector<double> psd_result(cpsd_result.size());

    for (size_t k = 0; k < cpsd_result.size(); ++k)
    {
        psd_result[k] = std::real(cpsd_result[k]);
    }

    return psd_result;
}

} // namespace msl::signal

#endif // MSL_POWER_SPECTRAL_DENSITY_HPP