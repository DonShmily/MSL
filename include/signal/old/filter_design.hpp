/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025
**  File: filter_design.hpp
**  Description: Digital filter design (IIR filters)
*/

#ifndef MSL_FILTER_DESIGN_HPP
#define MSL_FILTER_DESIGN_HPP

#include <cmath>
#include <complex>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>


namespace msl::filter
{

// ============================================================================
// Filter specifications
// ============================================================================

/**
 * @brief Filter type based on frequency response
 */
enum class FilterType
{
    lowpass,  // Low-pass filter
    highpass, // High-pass filter
    bandpass, // Band-pass filter
    bandstop  // Band-stop (notch) filter
};

/**
 * @brief Filter design method
 */
enum class FilterMethod
{
    butterworth, // Maximally flat passband
    chebyshev1,  // Equiripple passband, monotonic stopband (future)
    chebyshev2,  // Monotonic passband, equiripple stopband (future)
    elliptic     // Equiripple in both bands (future)
};

/**
 * @brief Filter coefficients (numerator and denominator)
 *
 * Represents digital filter as transfer function:
 * H(z) = (b[0] + b[1]*z^-1 + ... + b[M]*z^-M) / (a[0] + a[1]*z^-1 + ... +
 * a[N]*z^-N)
 *
 * Where a[0] is typically normalized to 1.0
 */
struct FilterCoefficients
{
    std::vector<double> b; // Numerator coefficients
    std::vector<double> a; // Denominator coefficients (a[0] = 1.0)

    FilterCoefficients() = default;
    FilterCoefficients(std::vector<double> num, std::vector<double> den)
        : b(std::move(num)), a(std::move(den))
    {}

    [[nodiscard]] size_t order() const
    {
        return std::max(a.size(), b.size()) - 1;
    }
};

// ============================================================================
// Butterworth filter design
// ============================================================================

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

        // Pre-warp cutoff frequency
        const double wc = std::tan(std::numbers::pi * fc);

        // Get analog prototype poles
        auto poles = get_analog_poles();

        // Apply frequency transformation and bilinear transform
        std::vector<std::complex<double>> digital_poles;
        digital_poles.reserve(order_);

        for (const auto &p : poles)
        {
            // Scale by cutoff frequency
            auto ps = p * wc;

            // Bilinear transformation: s -> (z-1)/(z+1)
            // pole in z-domain: (1 + ps)/(1 - ps)
            auto pz = (1.0 + ps) / (1.0 - ps);
            digital_poles.push_back(pz);
        }

        // Convert poles to polynomial coefficients
        coeffs_.a = poles_to_polynomial(digital_poles);

        // Compute numerator based on filter type
        if (type_ == FilterType::lowpass)
        {
            coeffs_.b = compute_lowpass_numerator();
        }
        else
        {
            coeffs_.b = compute_highpass_numerator();
        }

        // Normalize gain
        normalize_gain(fc);
    }

    /**
     * @brief Design bandpass filter
     */
    void design_bandpass()
    {
        // Center frequency and bandwidth
        const double w1 = std::tan(std::numbers::pi * fc_low_);
        const double w2 = std::tan(std::numbers::pi * fc_high_);
        const double w0 = std::sqrt(w1 * w2); // Geometric mean
        const double bw = w2 - w1;

        // Get analog prototype poles
        auto lp_poles = get_analog_poles();

        // Transform lowpass to bandpass (doubles the order)
        std::vector<std::complex<double>> bp_poles;
        bp_poles.reserve(2 * order_);

        for (const auto &p : lp_poles)
        {
            // Lowpass to bandpass transformation
            // s -> (s^2 + w0^2) / (s * bw)
            // Poles: s_bp = bw/2 * p ± sqrt((bw/2*p)^2 - w0^2)

            auto alpha = bw * p / 2.0;
            auto beta = std::sqrt(alpha * alpha - w0 * w0);

            bp_poles.push_back(alpha + beta);
            bp_poles.push_back(alpha - beta);
        }

        // Apply bilinear transformation
        std::vector<std::complex<double>> digital_poles;
        digital_poles.reserve(bp_poles.size());

        for (const auto &p : bp_poles)
        {
            auto pz = (1.0 + p) / (1.0 - p);
            digital_poles.push_back(pz);
        }

        // Convert to polynomial
        coeffs_.a = poles_to_polynomial(digital_poles);

        // Numerator for bandpass
        coeffs_.b.resize(2 * order_ + 1, 0.0);

        // Bandpass: numerator is [1, 0, -1, 0, 1, ...]
        for (int i = 0; i <= 2 * order_; i += 2)
        {
            coeffs_.b[i] = (i / 2) % 2 == 0 ? 1.0 : -1.0;
        }

        // Normalize gain at center frequency
        const double fc_center = (fc_low_ + fc_high_) / 2.0;
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
        std::vector<std::complex<double>> poly(poles.size() + 1, 0.0);
        poly[0] = 1.0;

        for (const auto &pole : poles)
        {
            // Multiply polynomial by (z - pole)
            for (size_t i = poly.size() - 1; i > 0; --i)
            {
                poly[i] = poly[i - 1] - pole * poly[i];
            }
            poly[0] *= -pole;
        }

        // Convert to real (imaginary parts should be ~0 for real filters)
        std::vector<double> result(poly.size());
        for (size_t i = 0; i < poly.size(); ++i)
        {
            result[i] = poly[i].real();
        }

        // Normalize so a[0] = 1
        double a0 = result[0];
        for (auto &coeff : result)
        {
            coeff /= a0;
        }

        return result;
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
        // Evaluate frequency response at frequency f
        const std::complex<double> j(0.0, 1.0);
        const std::complex<double> z = std::exp(j * 2.0 * std::numbers::pi * f);

        std::complex<double> num(0.0, 0.0);
        std::complex<double> den(0.0, 0.0);

        std::complex<double> zpower(1.0, 0.0);
        for (size_t i = 0; i < coeffs_.b.size(); ++i)
        {
            num += coeffs_.b[i] * zpower;
            zpower *= z;
        }

        zpower = std::complex<double>(1.0, 0.0);
        for (size_t i = 0; i < coeffs_.a.size(); ++i)
        {
            den += coeffs_.a[i] * zpower;
            zpower *= z;
        }

        double gain = std::abs(num / den);

        // Normalize numerator
        for (auto &coeff : coeffs_.b)
        {
            coeff /= gain;
        }
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

} // namespace msl::filter

#endif // MSL_FILTER_DESIGN_HPP