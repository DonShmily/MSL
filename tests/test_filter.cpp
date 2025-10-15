/**
 * Examples and tests for MSL filter module
 */

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>

#include "fft/fft.hpp" // For frequency analysis
#include "signal/butterworth_filter.hpp"
#include "signal/filter_apply.hpp"
#include "signal/filter_design.hpp"


#include "helper.hpp"

using namespace msl;
using namespace msl::signal;

// ============================================================================
// Test 1: Filter design - Butterworth lowpass
// ============================================================================

void test_butterworth_design()
{
    std::cout << "========================================\n";
    std::cout << "Test 1: Butterworth Filter Design\n";
    std::cout << "========================================\n\n";

    // Design 4th order lowpass filter
    // Cutoff at 0.2 (normalized frequency, 0.2 * Nyquist)
    ButterworthFilter filter(4, 0.2, FilterType::lowpass);

    print_coefficients(filter.coefficients(), "4th Order Lowpass (fc=0.2)");

    // Design highpass
    auto hp_coeffs = butterworth_highpass(4, 0.3);
    print_coefficients(hp_coeffs, "4th Order Highpass (fc=0.3)");

    // Design bandpass
    auto bp_coeffs = butterworth_bandpass(2, 0.1, 0.4);
    print_coefficients(bp_coeffs, "2nd Order Bandpass (0.1-0.4)");
}

// ============================================================================
// Test 2: Lowpass filtering
// ============================================================================

void test_lowpass_filter()
{
    std::cout << "========================================\n";
    std::cout << "Test 2: Lowpass Filtering\n";
    std::cout << "========================================\n\n";

    // Generate noisy signal: low freq + high freq noise
    const size_t n = 1000;
    const double fs = 1000.0; // 1 kHz sampling
    std::vector<double> signal(n);

    for (size_t i = 0; i < n; ++i)
    {
        double t = i / fs;
        // Signal: 10 Hz sine + 100 Hz noise
        signal[i] = std::sin(2.0 * std::numbers::pi * 10.0 * t)
                    + 0.5 * std::sin(2.0 * std::numbers::pi * 100.0 * t);
    }

    print_vector(signal, "Original signal (first 10)", 10);

    // Design filter: cutoff at 30 Hz
    // Normalized frequency: 30 / (fs/2) = 30 / 500 = 0.06
    double fc_norm = 30.0 / (fs / 2.0);

    // Apply zero-phase lowpass filter
    auto filtered = lowpass(signal, 4, fc_norm, true);
    print_vector(filtered, "Filtered signal (first 10)", 10);

    // Compute RMS before and after filtering
    auto compute_rms = [](const std::vector<double> &v) {
        double sum = 0.0;
        for (auto val : v)
            sum += val * val;
        return std::sqrt(sum / v.size());
    };

    std::cout << "\nRMS before: " << compute_rms(signal) << "\n";
    std::cout << "RMS after:  " << compute_rms(filtered) << "\n";
    std::cout << "Noise reduction: "
              << (1.0 - compute_rms(filtered) / compute_rms(signal)) * 100.0
              << "%\n\n";
}

// ============================================================================
// Test 3: Highpass filtering
// ============================================================================

void test_highpass_filter()
{
    std::cout << "========================================\n";
    std::cout << "Test 3: Highpass Filtering (DC Removal)\n";
    std::cout << "========================================\n\n";

    // Generate signal with DC offset
    const size_t n = 500;
    const double fs = 100.0;
    std::vector<double> signal(n);

    const double dc_offset = 5.0;
    for (size_t i = 0; i < n; ++i)
    {
        double t = i / fs;
        signal[i] = dc_offset + std::sin(2.0 * std::numbers::pi * 2.0 * t);
    }

    // Mean before filtering
    double mean_before = 0.0;
    for (auto val : signal)
        mean_before += val;
    mean_before /= signal.size();

    std::cout << "Signal mean before: " << mean_before << "\n";

    // Highpass filter: cutoff at 1 Hz (normalized: 1 / 50 = 0.02)
    auto filtered = highpass(signal, 4, 0.02, true);

    // Mean after filtering
    double mean_after = 0.0;
    for (auto val : filtered)
        mean_after += val;
    mean_after /= filtered.size();

    std::cout << "Signal mean after:  " << mean_after << "\n";
    std::cout << "DC removal: " << (dc_offset - mean_after) << "\n\n";
}

// ============================================================================
// Test 4: Bandpass filtering
// ============================================================================

void test_bandpass_filter()
{
    std::cout << "========================================\n";
    std::cout << "Test 4: Bandpass Filtering\n";
    std::cout << "========================================\n\n";

    // Generate multi-frequency signal
    const size_t n = 1024;
    const double fs = 1000.0;
    std::vector<double> signal(n);

    for (size_t i = 0; i < n; ++i)
    {
        double t = i / fs;
        // Three frequencies: 10 Hz, 50 Hz, 200 Hz
        signal[i] = std::sin(2.0 * std::numbers::pi * 10.0 * t)
                    + std::sin(2.0 * std::numbers::pi * 50.0 * t)
                    + std::sin(2.0 * std::numbers::pi * 200.0 * t);
    }

    std::cout << "Original signal contains: 10 Hz, 50 Hz, 200 Hz\n";

    // Bandpass filter: 30-100 Hz (normalized: 0.06-0.2)
    auto filtered = bandpass(signal, 4, 0.06, 0.2, true);

    // Analyze frequency content using FFT
    auto analyze = [&](const std::vector<double> &sig,
                       const std::string &name) {
        auto fft_result = fft::fft(sig);
        auto magnitude = fft::magnitude_spectrum(fft_result);
        auto freqs = fft::fft_frequencies(sig.size(), fs);

        // Find peaks
        std::cout << name << " frequency peaks:\n";
        for (size_t i = 1; i < magnitude.size() - 1; ++i)
        {
            if (magnitude[i] > 100.0 && // Threshold
                magnitude[i] > magnitude[i - 1]
                && magnitude[i] > magnitude[i + 1])
            {
                std::cout << "  " << freqs[i] << " Hz (mag: " << magnitude[i]
                          << ")\n";
            }
        }
    };

    analyze(signal, "Original");
    analyze(filtered, "Filtered");

    std::cout << "\nExpected: only 50 Hz component remains\n\n";
}

// ============================================================================
// Test 5: Forward vs Forward-Backward filtering
// ============================================================================

void test_phase_distortion()
{
    std::cout << "========================================\n";
    std::cout << "Test 5: Phase Distortion (filtfilt)\n";
    std::cout << "========================================\n\n";

    // Generate step function
    const size_t n = 200;
    std::vector<double> signal(n, 0.0);
    for (size_t i = n / 2; i < n; ++i)
    {
        signal[i] = 1.0;
    }

    // Design lowpass filter
    auto coeffs = butterworth_lowpass(4, 0.1);

    // Forward filtering only (has phase lag)
    auto forward = filter(signal, coeffs);

    // Forward-backward filtering (zero phase)
    auto zero_phase = filtfilt(signal, coeffs);

    // Find rise time (10% to 90%)
    auto find_rise_index = [](const std::vector<double> &v, double threshold) {
        for (size_t i = 0; i < v.size(); ++i)
        {
            if (v[i] >= threshold)
                return i;
        }
        return v.size();
    };

    size_t forward_50 = find_rise_index(forward, 0.5);
    size_t zerophase_50 = find_rise_index(zero_phase, 0.5);

    std::cout << "Step input at sample " << n / 2 << "\n";
    std::cout << "Forward filter reaches 50% at sample " << forward_50
              << " (delay: " << (int)forward_50 - (int)(n / 2) << ")\n";
    std::cout << "Zero-phase filter reaches 50% at sample " << zerophase_50
              << " (delay: " << (int)zerophase_50 - (int)(n / 2) << ")\n\n";
}

// ============================================================================
// Test 6: Matrix filtering (batch processing)
// ============================================================================

void test_matrix_filtering()
{
    std::cout << "========================================\n";
    std::cout << "Test 6: Matrix Filtering (Batch)\n";
    std::cout << "========================================\n\n";

    // Create matrix with multiple signals
    const size_t n_samples = 500;
    const size_t n_signals = 5;
    const double fs = 100.0;

    matrixd signals(n_samples, n_signals);

    // Each column is a signal at different frequency
    for (size_t j = 0; j < n_signals; ++j)
    {
        double freq = 5.0 * (j + 1); // 5, 10, 15, 20, 25 Hz
        for (size_t i = 0; i < n_samples; ++i)
        {
            double t = i / fs;
            signals(i, j) = std::sin(2.0 * std::numbers::pi * freq * t);
        }
    }

    std::cout << "Input: " << n_signals << " signals at 5, 10, 15, 20, 25 Hz\n";

    // Apply lowpass filter: cutoff at 12 Hz (0.24 normalized)
    auto coeffs = butterworth_lowpass(4, 0.24);
    auto filtered = filtfilt_columns(signals, coeffs);

    std::cout << "Applied lowpass filter (fc = 12 Hz)\n";
    std::cout << "Output matrix: " << filtered.rows() << "x" << filtered.cols()
              << "\n";

    // Check which signals pass
    for (size_t j = 0; j < n_signals; ++j)
    {
        double rms_in = 0.0, rms_out = 0.0;
        for (size_t i = 0; i < n_samples; ++i)
        {
            rms_in += signals(i, j) * signals(i, j);
            rms_out += filtered(i, j) * filtered(i, j);
        }
        rms_in = std::sqrt(rms_in / n_samples);
        rms_out = std::sqrt(rms_out / n_samples);

        double attenuation = rms_out / rms_in;
        std::cout << "Signal " << j << " (" << 5.0 * (j + 1) << " Hz): "
                  << "attenuation = " << std::setprecision(3) << attenuation
                  << (attenuation > 0.5 ? " (passed)" : " (blocked)") << "\n";
    }
    std::cout << "\n";
}

// ============================================================================
// Test 7: Frequency response
// ============================================================================

void test_frequency_response()
{
    std::cout << "========================================\n";
    std::cout << "Test 7: Filter Frequency Response\n";
    std::cout << "========================================\n\n";

    // Design filter
    auto coeffs = butterworth_lowpass(4, 0.2);

    std::cout << "Testing frequency response at various frequencies:\n";
    std::cout << std::setw(15) << "Freq (norm)" << std::setw(15) << "Magnitude"
              << std::setw(15) << "Magnitude (dB)\n";
    std::cout << std::string(45, '-') << "\n";

    for (double f = 0.05; f <= 0.5; f += 0.05)
    {
        // Evaluate H(e^jω)
        std::complex<double> j(0.0, 1.0);
        std::complex<double> z = std::exp(j * 2.0 * std::numbers::pi * f);

        std::complex<double> num(0.0, 0.0);
        std::complex<double> den(0.0, 0.0);

        std::complex<double> zpower(1.0, 0.0);
        for (size_t i = 0; i < coeffs.b.size(); ++i)
        {
            num += coeffs.b[i] * zpower;
            zpower *= z;
        }

        zpower = std::complex<double>(1.0, 0.0);
        for (size_t i = 0; i < coeffs.a.size(); ++i)
        {
            den += coeffs.a[i] * zpower;
            zpower *= z;
        }

        double mag = std::abs(num / den);
        double mag_db = 20.0 * std::log10(mag);

        std::cout << std::setw(15) << std::fixed << std::setprecision(3) << f
                  << std::setw(15) << std::setprecision(4) << mag
                  << std::setw(15) << std::setprecision(2) << mag_db << "\n";
    }
    std::cout << "\n";
}

// ============================================================================
// Test 8: Real-world example - ECG filtering
// ============================================================================

void test_ecg_filtering()
{
    std::cout << "========================================\n";
    std::cout << "Test 8: ECG Signal Filtering\n";
    std::cout << "========================================\n\n";

    // Simulate ECG signal (simplified)
    const size_t n = 1000;
    const double fs = 360.0; // 360 Hz (typical ECG sampling)
    std::vector<double> ecg(n);

    // ECG waveform (simplified)
    for (size_t i = 0; i < n; ++i)
    {
        double t = i / fs;
        // QRS complex (simplified as sharp peak)
        double qrs = 0.0;
        double t_mod = std::fmod(t, 0.8); // 75 bpm
        if (t_mod < 0.1)
        {
            qrs = 10.0 * std::exp(-100.0 * (t_mod - 0.05) * (t_mod - 0.05));
        }

        // Baseline wander (0.5 Hz)
        double baseline = 0.3 * std::sin(2.0 * std::numbers::pi * 0.5 * t);

        // 60 Hz power line noise
        double noise = 0.2 * std::sin(2.0 * std::numbers::pi * 60.0 * t);

        ecg[i] = qrs + baseline + noise;
    }

    std::cout << "Simulated ECG with baseline wander and 60 Hz noise\n\n";

    // Step 1: Remove baseline wander (highpass at 0.5 Hz)
    // Normalized: 0.5 / 180 = 0.0028
    auto step1 = highpass(ecg, 2, 0.0028, true);
    std::cout << "Step 1: Baseline wander removed (highpass 0.5 Hz)\n";

    // Step 2: Remove 60 Hz noise (lowpass at 40 Hz)
    // Normalized: 40 / 180 = 0.222
    auto step2 = lowpass(step1, 4, 0.222, true);
    std::cout << "Step 2: 60 Hz noise removed (lowpass 40 Hz)\n";

    // Alternatively: use bandpass directly (0.5 - 40 Hz)
    auto direct = bandpass(ecg, 4, 0.0028, 0.222, true);
    std::cout << "Direct: Bandpass filter (0.5 - 40 Hz)\n\n";

    // Compare SNR
    auto compute_snr = [&](const std::vector<double> &filtered) {
        // Simple SNR estimate
        double signal_power = 0.0, noise_power = 0.0;
        for (size_t i = 0; i < n; ++i)
        {
            signal_power += filtered[i] * filtered[i];
            noise_power += (ecg[i] - filtered[i]) * (ecg[i] - filtered[i]);
        }
        return 10.0 * std::log10(signal_power / noise_power);
    };

    std::cout << "SNR improvement (2-step): " << compute_snr(step2) << " dB\n";
    std::cout << "SNR improvement (direct): " << compute_snr(direct)
              << " dB\n\n";
}

// ============================================================================
// Main
// ============================================================================

int test_filter()
{
    try
    {
        test_butterworth_design();
        test_lowpass_filter();
        test_highpass_filter();
        test_bandpass_filter();
        test_phase_distortion();
        test_matrix_filtering();
        test_frequency_response();
        test_ecg_filtering();

        std::cout << "========================================\n";
        std::cout << "All filter tests completed successfully!\n";
        std::cout << "========================================\n";
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}