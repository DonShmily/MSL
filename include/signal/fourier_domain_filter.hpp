/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \include\signal\fourier_domain_filter.hpp
** -----
** File Created: Saturday, 18th October 2025 16:30:55
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Saturday, 18th October 2025 16:30:57
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_FOURIER_DOMAIN_FILTER_HPP
#define MSL_FOURIER_DOMAIN_FILTER_HPP

#include <cstddef>
#include <stdexcept>

#include "filter_design.hpp"

namespace msl::signal
{
/**
 * @brief Fourier domain filter designer
 *
 * Designs bandpass or bandstop filters in the Fourier domain.
 */
class FourierDomainFilter
{
private:
    double fc_low_{0.0};  // Low cutoff frequency (0-1 where 1 = Nyquist)
    double fc_high_{1.0}; // High cutoff frequency (0-1 where 1 = Nyquist)
    std::size_t nfft_{0}; // FFT size, 0 = auto

    FilterType type_;

public:
    FourierDomainFilter() = default;
    /**
     * @brief Construct Fourier domain filter
     *
     * @param fc_low Low cutoff frequency (normalized: 0 < fc < 1, where 1 =
     * Nyquist frequency)
     * @param fc_high High cutoff frequency (normalized: 0 < fc < 1, where 1 =
     * Nyquist frequency)
     * @param nfft FFT size
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
        if (fc_low <= 0.0 || fc_low >= 1.0 || fc_high <= 0.0 || fc_high >= 1.0)
        {
            throw std::invalid_argument("Frequencies must be normalized (0 < f "
                                        "< 1, where 1 = Nyquist)");
        }

        if (fc_low >= fc_high)
        {
            throw std::invalid_argument(
                "Low frequency must be < high frequency");
        }
    }
};

} // namespace msl::signal

#endif // MSL_FOURIER_DOMAIN_FILTER_HPP