/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: \include\signal\butterworth_filter.hpp
** -----
** File Created: Wednesday, 15th October 2025 14:16:08
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Wednesday, 15th October 2025 14:16:20
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/
#ifndef MSL_BUTTERWORTH_FILTER_HPP
#define MSL_BUTTERWORTH_FILTER_HPP

#include "filter_apply.hpp"
#include "filter_design.hpp"

#include <cmath>
#include <complex>

namespace msl::signal
{
/**
 * @brief Butterworth filter designer
 *
 * Butterworth filters have maximally flat magnitude response in the passband.
 * No ripple in passband or stopband, monotonic response.
 *
 * The digital filter is designed using bilinear transformation from analog
 * prototype.
 */
class ButterworthFilter
{
private:
    int order_;
    double fc_low_;  // Low cutoff frequency (normalized, 0-1 where 1 = Nyquist)
    double fc_high_; // High cutoff frequency
    FilterType type_;

    FilterCoefficients coeffs_;

public:
    /**
     * @brief Default constructor
     */
    ButterworthFilter() = default;

    /**
     * @brief Construct lowpass/highpass filter
     *
     * @param order Filter order (must be positive)
     * @param fc Cutoff frequency (normalized: 0 < fc < 1, where 1 = Nyquist
     * frequency)
     * @param type Filter type (lowpass or highpass)
     *
     * @throws std::invalid_argument if parameters are invalid
     */
    ButterworthFilter(int order,
                      double fc,
                      FilterType type = FilterType::lowpass)
        : order_(order), type_(type)
    {
        validate_order(order);
        validate_frequency(fc);

        if (type == FilterType::lowpass)
        {
            fc_low_ = fc;
            fc_high_ = 0.0;
        }
        else if (type == FilterType::highpass)
        {
            fc_low_ = 0.0;
            fc_high_ = fc;
        }
        else
        {
            throw std::invalid_argument(
                "Single frequency constructor only for lowpass/highpass");
        }

        design();
    }

    /**
     * @brief Construct bandpass/bandstop filter
     *
     * @param order Filter order (must be positive)
     * @param fc_low Low cutoff frequency (normalized)
     * @param fc_high High cutoff frequency (normalized)
     * @param type Filter type (bandpass or bandstop)
     *
     * @throws std::invalid_argument if parameters are invalid
     */
    ButterworthFilter(int order,
                      double fc_low,
                      double fc_high,
                      FilterType type = FilterType::bandpass)
        : order_(order), fc_low_(fc_low), fc_high_(fc_high), type_(type)
    {
        validate_order(order);
        validate_frequency(fc_low);
        validate_frequency(fc_high);

        if (fc_low >= fc_high)
        {
            throw std::invalid_argument(
                "Low frequency must be < high frequency");
        }

        design();
    }

    /**
     * @brief Get filter coefficients
     */
    [[nodiscard]] const FilterCoefficients &coefficients() const
    {
        return coeffs_;
    }

    /**
     * @brief Get numerator coefficients
     */
    [[nodiscard]] const std::vector<double> &b() const { return coeffs_.b; }

    /**
     * @brief Get denominator coefficients
     */
    [[nodiscard]] const std::vector<double> &a() const { return coeffs_.a; }

    /**
     * @brief Get filter specifications
     */
    [[nodiscard]] int order() const { return order_; }
    [[nodiscard]] FilterType type() const { return type_; }
    [[nodiscard]] double fc_low() const { return fc_low_; }
    [[nodiscard]] double fc_high() const { return fc_high_; }

    /**
     * @brief Redesign filter with new parameters
     */
    void redesign(int order, double fc, FilterType type = FilterType::lowpass)
    {
        validate_order(order);
        validate_frequency(fc);

        order_ = order;
        type_ = type;

        if (type == FilterType::lowpass)
        {
            fc_low_ = fc;
            fc_high_ = 0.0;
        }
        else
        {
            fc_low_ = 0.0;
            fc_high_ = fc;
        }

        design();
    }

    void redesign(int order,
                  double fc_low,
                  double fc_high,
                  FilterType type = FilterType::bandpass)
    {
        validate_order(order);
        validate_frequency(fc_low);
        validate_frequency(fc_high);

        if (fc_low >= fc_high)
        {
            throw std::invalid_argument(
                "Low frequency must be < high frequency");
        }

        order_ = order;
        fc_low_ = fc_low;
        fc_high_ = fc_high;
        type_ = type;

        design();
    }

private:
    void validate_order(int order) const
    {
        if (order <= 0)
        {
            throw std::invalid_argument("Filter order must be positive");
        }
    }

    void validate_frequency(double f) const
    {
        if (f <= 0.0 || f >= 1.0)
        {
            throw std::invalid_argument(
                "Frequency must be normalized (0 < f < 1, where 1 = Nyquist)");
        }
    }

    /**
     * @brief Design the filter
     */
    void design()
    {
        if (type_ == FilterType::bandpass || type_ == FilterType::bandstop)
        {
            design_bandpass();
        }
        else
        {
            design_lowpass_highpass();
        }
    }

    /**
     * @brief Design lowpass or highpass filter
     */
    void design_lowpass_highpass()
    {
        const double fc = (type_ == FilterType::lowpass) ? fc_low_ : fc_high_;

        // Pre-warp cutoff frequency (T = 1 assumed)
        const double wc =
            2.0 * std::tan(std::numbers::pi * fc / 2.0); // <-- 修正

        // Analog prototype poles (normalized)
        auto analog_poles = get_analog_poles();

        // Map analog poles to digital poles using bilinear transform:
        // s = p * wc, z = (2 + s) / (2 - s)
        std::vector<std::complex<double>> digital_poles;
        digital_poles.reserve(order_);

        for (const auto &p : analog_poles)
        {
            std::complex<double> s = p * wc; // analog pole scaled
            std::complex<double> pz =
                (2.0 + s) / (2.0 - s); // <-- 修正: T=1 mapping
            digital_poles.push_back(pz);
        }

        // Denominator polynomial in z^{-1}: prod (1 - pz * z^{-1})
        coeffs_.a =
            poles_to_polynomial(digital_poles); // now returns a in z^{-k} order

        // Numerator: for lowpass the digital numerator corresponds to (1 +
        // z^{-1})^N (binomial)
        if (type_ == FilterType::lowpass)
        {
            coeffs_.b = compute_lowpass_numerator(); // returns coefficients in
                                                     // z^{-k} order (C(N,k))
        }
        else // highpass
        {
            coeffs_.b = compute_highpass_numerator(); // also z^{-k} order
        }

        // Normalize: choose frequency point according to type
        if (type_ == FilterType::lowpass)
            normalize_gain(0.0); // unity at DC
        else
            normalize_gain(1.0); // unity at Nyquist (f=1 corresponds to ω=π)
    }

    /**
     * @brief Design bandpass filter
     */
    void design_bandpass()
    {
        // Pre-warp digital cutoff freqs (T = 1)
        const double wc1 = 2.0 * std::tan(std::numbers::pi * fc_low_ / 2.0);
        const double wc2 = 2.0 * std::tan(std::numbers::pi * fc_high_ / 2.0);

        const double w0 = std::sqrt(wc1 * wc2); // analog center freq
        const double bw = wc2 - wc1;            // analog bandwidth

        // Get analog lowpass prototype poles (order_ of them)
        auto lp_poles = get_analog_poles(); // returns poles of normalized LP

        // Transform each lowpass pole p -> two bandpass poles
        std::vector<std::complex<double>> bp_poles;
        bp_poles.reserve(2 * order_);

        for (const auto &p : lp_poles)
        {
            // alpha = (bw/2) * p
            std::complex<double> alpha = (bw / 2.0) * p;
            // beta = sqrt(alpha^2 - w0^2)
            std::complex<double> beta =
                std::sqrt(alpha * alpha - std::complex<double>(w0 * w0, 0.0));
            // two bandpass analog poles
            std::complex<double> s1 = alpha + beta;
            std::complex<double> s2 = alpha - beta;
            bp_poles.push_back(s1);
            bp_poles.push_back(s2);
        }

        // Bilinear transform each analog s pole to z-domain: z = (2 + s) / (2 -
        // s)
        std::vector<std::complex<double>> digital_poles;
        digital_poles.reserve(bp_poles.size());
        for (const auto &s : bp_poles)
        {
            std::complex<double> z = (2.0 + s) / (2.0 - s);
            digital_poles.push_back(z);
        }

        // Denominator polynomial from digital poles (z^{-1} representation)
        coeffs_.a = poles_to_polynomial(
            digital_poles); // returns [a0,a1,...,aM] for z^{-1} powers

        // Numerator: zeros of analog BP at s=0 (multiplicity N) map to z=1
        // (multiplicity N), and zeros at s=infty map to z=-1 (multiplicity N)
        // -> so digital numerator has factor (1 - z^{-2})^N. Build b of length
        // 2*N + 1 with only even coefficients:
        const int M = 2 * order_;
        coeffs_.b.assign(M + 1, 0.0);

        // binomial coefficients for (1 - t^2)^N where t = z^{-1}
        // (1 - t^2)^N = sum_{m=0..N} C(N,m) * (-1)^m * t^{2m}
        auto binom = [](int n, int k) -> double {
            double r = 1.0;
            for (int i = 1; i <= k; ++i)
                r *= double(n - (k - i)) / double(i);
            return r;
        };

        for (int m = 0; m <= order_; ++m)
        {
            double c = binom(order_, m)
                       * ((m % 2 == 0) ? 1.0 : -1.0); // (-1)^m * C(N,m)
            coeffs_.b[2 * m] = c;
        }

        // Normalize gain at center frequency (use normalized digital freq
        // f_center in 0..1)
        const double fc_center =
            0.5 * (fc_low_ + fc_high_); // normalized (0..1)
        normalize_gain(fc_center);
    }


    /**
     * @brief Get analog Butterworth poles
     *
     * Returns poles of normalized lowpass Butterworth filter (cutoff = 1 rad/s)
     */
    std::vector<std::complex<double>> get_analog_poles() const
    {
        std::vector<std::complex<double>> poles;
        poles.reserve(order_);

        for (int k = 0; k < order_; ++k)
        {
            // Pole angle: (2k + 1)π / (2N) + π/2
            double angle = std::numbers::pi * (2.0 * k + 1.0) / (2.0 * order_)
                           + std::numbers::pi / 2.0;

            poles.emplace_back(std::cos(angle), std::sin(angle));
        }

        return poles;
    }

    /**
     * @brief Convert poles to polynomial coefficients
     *
     * Given poles p1, p2, ..., pN, compute coefficients of:
     * (z - p1)(z - p2)...(z - pN) = a[0] + a[1]*z + ... + a[N]*z^N
     */
    std::vector<double>
    poles_to_polynomial(const std::vector<std::complex<double>> &poles) const
    {
        // We compute polynomial in z^{-1} form:
        //   A(z) = ∏ (1 - p_i z^{-1}) = a0 + a1 z^{-1} + ... + aN z^{-N}
        // Start with poly = [1]
        std::vector<std::complex<double>> poly(1);
        poly[0] = std::complex<double>(1.0, 0.0);

        for (const auto &p : poles)
        {
            std::vector<std::complex<double>> next(
                poly.size() + 1, std::complex<double>(0.0, 0.0));
            for (size_t k = 0; k < poly.size(); ++k)
            {
                // multiply existing terms by 1 (coefficient for z^0) -> shift 0
                next[k] += poly[k];
                // multiply existing terms by (-p) and shift by one power
                // (z^{-1})
                next[k + 1] += -p * poly[k];
            }
            poly.swap(next);
        }

        // Convert to real (imag parts should be ~0)
        std::vector<double> result(poly.size());
        for (size_t i = 0; i < poly.size(); ++i)
            result[i] = poly[i].real();

        // Normalize so a[0] == 1.0 (should already be 1)
        double a0 = result[0];
        if (std::abs(a0) < 1e-300)
            throw std::runtime_error(
                "numerical instability in poles_to_polynomial");
        for (auto &c : result)
            c /= a0;

        return result; // coefficients [a0, a1, ..., aN] corresponding to
                       // z^{-0}, z^{-1}, ...
    }

    /**
     * @brief Compute numerator for lowpass filter
     */
    std::vector<double> compute_lowpass_numerator() const
    {
        std::vector<double> b(order_ + 1);

        // Binomial coefficients
        b[0] = 1.0;
        for (int i = 1; i <= order_; ++i)
        {
            b[i] = b[i - 1] * (order_ - i + 1) / i;
        }

        return b;
    }

    /**
     * @brief Compute numerator for highpass filter
     */
    std::vector<double> compute_highpass_numerator() const
    {
        auto b = compute_lowpass_numerator();

        // Alternate signs for highpass
        for (int i = 0; i <= order_; ++i)
        {
            if (i % 2 == 1)
            {
                b[i] = -b[i];
            }
        }

        return b;
    }

    /**
     * @brief Normalize filter gain at specified frequency
     */
    void normalize_gain(double f)
    {
        // Evaluate H(e^{jω}) at frequency f where f in [0,1], 1 => Nyquist => ω
        // = π
        const std::complex<double> j(0.0, 1.0);
        const double omega = std::numbers::pi * f; // ω = π * f
        const std::complex<double> inv_z =
            std::exp(-j * omega); // z^{-1} = e^{-j ω}

        // Evaluate numerator and denominator at z^{-k}
        std::complex<double> num(0.0, 0.0), den(0.0, 0.0);

        std::complex<double> zpow(1.0, 0.0); // z^0 = 1
        for (size_t k = 0; k < coeffs_.b.size(); ++k)
        {
            num += coeffs_.b[k] * zpow;
            zpow *= inv_z;
        }

        zpow = std::complex<double>(1.0, 0.0);
        for (size_t k = 0; k < coeffs_.a.size(); ++k)
        {
            den += coeffs_.a[k] * zpow;
            zpow *= inv_z;
        }

        double gain = std::abs(num / den);
        if (gain == 0.0)
            return; // avoid division by zero

        for (auto &c : coeffs_.b)
            c /= gain;
    }
};

// ============================================================================
// Convenience functions
// ============================================================================

/**
 * @brief Design lowpass Butterworth filter
 *
 * @param order Filter order
 * @param fc Cutoff frequency (normalized, 0 < fc < 1)
 * @return Filter coefficients
 */
inline FilterCoefficients butterworth_lowpass(int order, double fc)
{
    return ButterworthFilter(order, fc, FilterType::lowpass).coefficients();
}

/**
 * @brief Design highpass Butterworth filter
 */
inline FilterCoefficients butterworth_highpass(int order, double fc)
{
    return ButterworthFilter(order, fc, FilterType::highpass).coefficients();
}

/**
 * @brief Design bandpass Butterworth filter
 */
inline FilterCoefficients
butterworth_bandpass(int order, double fc_low, double fc_high)
{
    return ButterworthFilter(order, fc_low, fc_high, FilterType::bandpass)
        .coefficients();
}

// ============================================================================
// Convenience functions
// ============================================================================

/**
 * @brief Design and apply lowpass filter in one step
 *
 * @param signal Input signal
 * @param order Filter order
 * @param fc Cutoff frequency (normalized)
 * @param zero_phase Use zero-phase filtering (default: true)
 * @return Filtered signal
 */
inline std::vector<double> lowpass(std::span<const double> signal,
                                   int order,
                                   double fc,
                                   bool zero_phase = true)
{
    auto coeffs = butterworth_lowpass(order, fc);
    return zero_phase ? filtfilt(signal, coeffs) : filter(signal, coeffs);
}

/**
 * @brief Design and apply highpass filter
 */
inline std::vector<double> highpass(std::span<const double> signal,
                                    int order,
                                    double fc,
                                    bool zero_phase = true)
{
    auto coeffs = butterworth_highpass(order, fc);
    return zero_phase ? filtfilt(signal, coeffs) : filter(signal, coeffs);
}

/**
 * @brief Design and apply bandpass filter
 */
inline std::vector<double> bandpass(std::span<const double> signal,
                                    int order,
                                    double fc_low,
                                    double fc_high,
                                    bool zero_phase = true)
{
    auto coeffs = butterworth_bandpass(order, fc_low, fc_high);
    return zero_phase ? filtfilt(signal, coeffs) : filter(signal, coeffs);
}

// Vector overloads
inline std::vector<double> lowpass(const std::vector<double> &signal,
                                   int order,
                                   double fc,
                                   bool zero_phase = true)
{
    return lowpass(std::span<const double>(signal), order, fc, zero_phase);
}

inline std::vector<double> highpass(const std::vector<double> &signal,
                                    int order,
                                    double fc,
                                    bool zero_phase = true)
{
    return highpass(std::span<const double>(signal), order, fc, zero_phase);
}

inline std::vector<double> bandpass(const std::vector<double> &signal,
                                    int order,
                                    double fc_low,
                                    double fc_high,
                                    bool zero_phase = true)
{
    return bandpass(
        std::span<const double>(signal), order, fc_low, fc_high, zero_phase);
}

} // namespace msl::signal

#endif // MSL_BUTTERWORTH_FILTER_HPP