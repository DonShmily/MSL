/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \include\signal\fourier_domain_filter.hpp
** -----
** File Created: Saturday, 18th October 2025
** Author: Dong Feiyue (FeiyueDong@outlook.com)
*/
#ifndef MSL_FOURIER_DOMAIN_FILTER_HPP
#define MSL_FOURIER_DOMAIN_FILTER_HPP

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <functional>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "fft.hpp"
#include "filter_design.hpp"

namespace msl::signal
{
/**
 * @brief Fourier domain filter designer
 *
 * Designs bandpass or bandstop filters in the Fourier domain.
 * The filter is designed by creating an ideal frequency response and
 * applying a window function to reduce ringing artifacts.
 */
class FourierDomainFilter
{
public:
    enum class WindowType
    {
        rectangular,
        hamming,
        hanning,
        blackman
    };

private:
    double fc_low_{0.0};  // Low cutoff frequency (0-1 where 1 = Nyquist)
    double fc_high_{1.0}; // High cutoff frequency (0-1 where 1 = Nyquist)
    std::size_t nfft_{0}; // FFT size, 0 = auto

    FilterType type_{FilterType::bandpass};

    WindowType window_type_{WindowType::rectangular};
    double transition_band_{0.0}; // Transition band width (normalized)

    std::function<double(double)> window_function_;

public:
    FourierDomainFilter() = default;

    /**
     * @brief Construct Fourier domain filter
     *
     * @param fc_low Low cutoff frequency (normalized: 0 < fc < 1, where 1 =
     * Nyquist frequency)
     * @param fc_high High cutoff frequency (normalized: 0 < fc < 1, where 1 =
     * Nyquist frequency)
     * @param nfft FFT size (0 = auto-determine based on transition band)
     * @param type Filter type (bandpass or bandstop)
     *
     * @throws std::invalid_argument if parameters are invalid
     */
    FourierDomainFilter(double fc_low,
                        double fc_high,
                        std::size_t nfft = 0,
                        FilterType type = FilterType::bandpass)
        : fc_low_(fc_low), fc_high_(fc_high), nfft_(nfft), type_(type)
    {
        validate_parameters();
        design();
    }

    /**
     * @brief Set window type
     */
    void set_window_type(WindowType window_type)
    {
        window_type_ = window_type;
        design();
    }

    /**
     * @brief Set transition band width
     * @param width Transition band width (normalized, e.g., 0.01 = 1% of
     * Nyquist)
     */
    void set_transition_band(double width)
    {
        if (width < 0.0 || width >= 0.5)
        {
            throw std::invalid_argument(
                "Transition band must be in range (0, 0.5)");
        }
        transition_band_ = width;
        design();
    }

    /**
     * @brief Redesign the filter with new parameters
     */
    void redesign(double fc_low,
                  double fc_high,
                  FilterType type = FilterType::bandpass)
    {
        fc_low_ = fc_low;
        fc_high_ = fc_high;
        type_ = type;
        validate_parameters();
        design();
    }

    /**
     * @brief Apply the filter to the input signal
     *
     * @param signal Input signal (time domain)
     * @return Filtered signal (time domain)
     */
    std::vector<double> apply(const std::vector<double> &signal) const
    {
        std::size_t fft_size = (nfft_ == 0) ? signal.size() : nfft_;
        auto fft_data = fft::fft(signal, fft_size);
        std::vector<std::complex<double>> filtered_fft(fft_data.size());

        for (std::size_t i = 0; i < fft_data.size(); ++i)
        {
            double freq = static_cast<double>(i) / fft_size;
            if (freq > 0.5)
                freq -= 1.0; // Map to [-0.5, 0.5]
            freq *= 2.0;     // Normalize to [-1, 1]
            double gain = window_function_(freq);
            if (type_ == FilterType::bandstop)
            {
                gain = 1.0 - gain;
            }
            filtered_fft[i] = fft_data[i] * gain;
        }

        return fft::ifft_real(filtered_fft);
    }


private:
    void validate_parameters() const
    {
        if (fc_low_ <= 0.0 || fc_low_ >= 1.0 || fc_high_ <= 0.0
            || fc_high_ >= 1.0)
        {
            throw std::invalid_argument("Frequencies must be normalized (0 < f "
                                        "< 1, where 1 = Nyquist)");
        }

        if (fc_low_ >= fc_high_)
        {
            throw std::invalid_argument(
                "Low frequency must be < high frequency");
        }
    }

    /**
     * @brief Design the filter
     */
    void design()
    {
        window_function_ = [=, this](double f) -> double {
            f = std::abs(f); // Use absolute frequency
            const double f_low_start =
                std::max(0.0, fc_low_ - transition_band_);
            const double f_low_end = std::min(1.0, fc_low_ + transition_band_);
            const double f_high_start =
                std::max(0.0, fc_high_ - transition_band_);
            const double f_high_end =
                std::min(1.0, fc_high_ + transition_band_);


            // 全阻带
            if (f < f_low_start || f > f_high_end)
                return 0.0;
            // 全通带
            if (f >= f_low_end && f <= f_high_start)
                return 1.0;

            // 过渡带计算
            auto smooth = [&](double x) {
                double t = std::max(
                    std::min((x + transition_band_) / (2 * transition_band_),
                             1.0),
                    0.0);
                switch (window_type_)
                {
                    case WindowType::rectangular:
                        return (t >= 0.5) ? 1.0 : 0.0;
                    case WindowType::hanning:
                        return 0.5 * (1 - std::cos(std::numbers::pi * t));
                    case WindowType::hamming:
                        return 0.54 - 0.46 * std::cos(std::numbers::pi * t);
                    case WindowType::blackman:
                        return 0.42 - 0.5 * std::cos(std::numbers::pi * t)
                               + 0.08 * std::cos(2 * std::numbers::pi * t);
                }
                return 0.0;
            };

            if (f >= f_low_start && f < f_low_end)
                return smooth(f - f_low_start);
            if (f > f_high_start && f <= f_high_end)
                return smooth(f_high_end - f);
            return 1.0;
        };
    }
};

// ============================================================================
// Convenience functions
// ============================================================================

/**
 * @brief Design and apply bandpass filter in Fourier domain
 *
 * @param signal Input signal
 * @param fc_low Low cutoff frequency (normalized, 0 < fc < 1)
 * @param fc_high High cutoff frequency (normalized, 0 < fc < 1)
 * @param window_type Window function type
 * @return Filtered signal
 */
inline std::vector<double>
fourier_bandpass(const std::vector<double> &signal,
                 double fc_low,
                 double fc_high,
                 FourierDomainFilter::WindowType window_type =
                     FourierDomainFilter::WindowType::rectangular)
{
    FourierDomainFilter filter(fc_low, fc_high, 0, FilterType::bandpass);
    filter.set_window_type(window_type);
    return filter.apply(signal);
}

/**
 * @brief Design and apply bandstop filter in Fourier domain
 *
 * @param signal Input signal
 * @param fc_low Low cutoff frequency (normalized, 0 < fc < 1)
 * @param fc_high High cutoff frequency (normalized, 0 < fc < 1)
 * @param window_type Window function type
 * @return Filtered signal
 */
inline std::vector<double>
fourier_bandstop(const std::vector<double> &signal,
                 double fc_low,
                 double fc_high,
                 FourierDomainFilter::WindowType window_type =
                     FourierDomainFilter::WindowType::rectangular)
{
    FourierDomainFilter filter(fc_low, fc_high, 0, FilterType::bandstop);
    filter.set_window_type(window_type);
    return filter.apply(signal);
}

} // namespace msl::signal

#endif // MSL_FOURIER_DOMAIN_FILTER_HPP