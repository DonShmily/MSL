#include "helper.hpp"
#include "test.hpp"

// ============================================================================
// Test 1: Butterworth Lowpass Filter Coefficients
// ============================================================================

void test1_butterworth_coefficients()
{
    std::cout << "\n========================================\n";
    std::cout << "Test 1: Butterworth Filter Coefficients\n";
    std::cout << "========================================\n\n";

    // Test various configurations
    struct TestConfig
    {
        std::string name;
        int order;
        double fc_norm;
        signal::FilterType type;
    };

    std::vector<TestConfig> configs = {
        {"2nd order lowpass (fc=0.2)", 2, 0.2, signal::FilterType::lowpass},
        {"4th order lowpass (fc=0.2)", 4, 0.2, signal::FilterType::lowpass},
        {"4th order highpass (fc=0.3)", 4, 0.3, signal::FilterType::highpass},
    };

    for (const auto &config : configs)
    {
        std::cout << "Config: " << config.name << "\n";

        signal::ButterworthFilter filt(
            config.order, config.fc_norm, config.type);
        auto coeffs = filt.coefficients();

        std::cout << "  b coefficients: [";
        for (size_t i = 0; i < coeffs.b.size(); ++i)
        {
            std::cout << std::setprecision(6) << coeffs.b[i];
            if (i < coeffs.b.size() - 1)
                std::cout << ", ";
        }
        std::cout << "]\n";

        std::cout << "  a coefficients: [";
        for (size_t i = 0; i < coeffs.a.size(); ++i)
        {
            std::cout << std::setprecision(6) << coeffs.a[i];
            if (i < coeffs.a.size() - 1)
                std::cout << ", ";
        }
        std::cout << "]\n\n";
    }

    // Try to load and compare MATLAB results
    try
    {
        auto matlab_b =
            read_csv("msl_validation_data/test5_filter_b_coeffs.csv");
        auto matlab_a =
            read_csv("msl_validation_data/test5_filter_a_coeffs.csv");

        // Design same filter as MATLAB test 5
        auto coeffs = signal::butterworth_lowpass(4, 0.06); // 30 Hz / 500 Hz

        double error_b = max_abs_error(coeffs.b, matlab_b);
        double error_a = max_abs_error(coeffs.a, matlab_a);

        print_result("Filter b coefficients", error_b, 1e-10);
        print_result("Filter a coefficients", error_a, 1e-10);
    }
    catch (const std::exception &e)
    {
        std::cout << "Could not load MATLAB reference: " << e.what() << "\n";
    }
}

// ============================================================================
// Test 2: Lowpass Filtering
// ============================================================================

void test2_lowpass_filtering()
{
    std::cout << "\n========================================\n";
    std::cout << "Test 2: Lowpass Filtering\n";
    std::cout << "========================================\n\n";

    // Parameters (match MATLAB test 5)
    const size_t n = 1000;
    const double fs = 1000.0;
    const double fc_hz = 30.0;
    const double fc_norm = fc_hz / (fs / 2.0);

    // Generate signal: 10 Hz + 100 Hz noise
    std::vector<double> signal(n);
    for (size_t i = 0; i < n; ++i)
    {
        double t = i / fs;
        signal[i] = std::sin(2.0 * std::numbers::pi * 10.0 * t)
                    + 0.5 * std::sin(2.0 * std::numbers::pi * 100.0 * t);
    }

    // Apply lowpass filter
    auto filtered = signal::lowpass(signal, 4, fc_norm, true);

    // Compute RMS before and after
    auto compute_rms = [](const std::vector<double> &v) {
        double sum = 0.0;
        for (auto val : v)
            sum += val * val;
        return std::sqrt(sum / v.size());
    };

    double rms_before = compute_rms(signal);
    double rms_after = compute_rms(filtered);

    std::cout << "RMS before filtering: " << rms_before << "\n";
    std::cout << "RMS after filtering:  " << rms_after << "\n";
    std::cout << "Noise reduction: " << (1.0 - rms_after / rms_before) * 100.0
              << "%\n";

    // Compare with MATLAB
    try
    {
        auto matlab_signal =
            read_csv("msl_validation_data/test5_original_signal.csv");
        auto matlab_filtered =
            read_csv("msl_validation_data/test5_filtered_zerophase.csv");

        double signal_error = max_abs_error(signal, matlab_signal);
        double filtered_error = max_abs_error(filtered, matlab_filtered);

        print_result("Input signal match", signal_error, 1e-10);
        print_result("Filtered signal match", filtered_error, 1e-6);
    }
    catch (const std::exception &e)
    {
        std::cout << "\nCould not load MATLAB reference: " << e.what() << "\n";
    }
}

// ============================================================================
// Test 3: Highpass Filtering (DC Removal)
// ============================================================================

void test3_highpass_filtering()
{
    std::cout << "\n========================================\n";
    std::cout << "Test 3: Highpass Filtering\n";
    std::cout << "========================================\n\n";

    const size_t n = 500;
    const double fs = 100.0;
    const double dc_offset = 5.0;

    // Generate signal with DC offset
    std::vector<double> signal(n);
    for (size_t i = 0; i < n; ++i)
    {
        double t = i / fs;
        signal[i] = dc_offset + std::sin(2.0 * std::numbers::pi * 2.0 * t);
    }

    // Compute mean before
    double mean_before = 0.0;
    for (auto val : signal)
        mean_before += val;
    mean_before /= signal.size();

    // Apply highpass filter (1 Hz cutoff)
    auto filtered = signal::highpass(signal, 4, 0.02, true); // 1 / 50 = 0.02

    // Compute mean after
    double mean_after = 0.0;
    for (auto val : filtered)
        mean_after += val;
    mean_after /= filtered.size();

    std::cout << "Signal mean before: " << mean_before << "\n";
    std::cout << "Signal mean after:  " << mean_after << "\n";
    std::cout << "DC removed: " << (dc_offset - mean_after) << "\n";

    if (std::abs(mean_after) < 0.1)
    {
        std::cout << "✓ DC component successfully removed\n";
    }
    else
    {
        std::cout << "✗ DC removal incomplete\n";
    }
}

// ============================================================================
// Test 4: Bandpass Filtering
// ============================================================================

void test4_bandpass_filtering()
{
    std::cout << "\n========================================\n";
    std::cout << "Test 4: Bandpass Filtering\n";
    std::cout << "========================================\n\n";

    const size_t n = 1024;
    const double fs = 1000.0;

    // Generate signal: 10, 50, 200 Hz
    std::vector<double> signal(n);
    for (size_t i = 0; i < n; ++i)
    {
        double t = i / fs;
        signal[i] = std::sin(2.0 * std::numbers::pi * 10.0 * t)
                    + std::sin(2.0 * std::numbers::pi * 50.0 * t)
                    + std::sin(2.0 * std::numbers::pi * 200.0 * t);
    }

    // Apply bandpass filter: 30-100 Hz
    auto filtered = signal::bandpass(signal, 4, 0.06, 0.2, true);

    // Analyze frequency content
    auto analyze = [&](const std::vector<double> &sig,
                       const std::string &name) {
        auto fft_result = fft::fft(sig);
        auto magnitude = fft::magnitude_spectrum(fft_result);
        auto freqs = fft::fft_frequencies(sig.size(), fs);

        std::cout << name << " frequency peaks:\n";
        for (size_t i = 1; i < magnitude.size() - 1; ++i)
        {
            if (magnitude[i] > 100.0 && magnitude[i] > magnitude[i - 1]
                && magnitude[i] > magnitude[i + 1])
            {
                std::cout << "  " << freqs[i] << " Hz (mag: " << magnitude[i]
                          << ")\n";
            }
        }
    };

    analyze(signal, "Original");
    analyze(filtered, "Filtered");

    std::cout << "\nExpected: Only 50 Hz component should remain\n";
}

// ============================================================================
// Test 5: Phase Distortion (filtfilt vs single-pass)
// ============================================================================

void test5_phase_distortion()
{
    std::cout << "\n========================================\n";
    std::cout << "Test 5: Phase Distortion\n";
    std::cout << "========================================\n\n";

    // Generate step function
    const size_t n = 200;
    std::vector<double> signal(n, 0.0);
    for (size_t i = n / 2; i < n; ++i)
    {
        signal[i] = 1.0;
    }

    // Design filter
    auto coeffs = signal::butterworth_lowpass(4, 0.1);

    // Single-pass filtering
    auto forward = signal::filter(signal, coeffs);

    // Zero-phase filtering
    auto zero_phase = signal::filtfilt(signal, coeffs);

    // Find 50% rise points
    auto find_50_percent = [](const std::vector<double> &v) -> size_t {
        for (size_t i = 0; i < v.size(); ++i)
        {
            if (v[i] >= 0.5)
                return i;
        }
        return v.size();
    };

    size_t step_idx = n / 2;
    size_t forward_50 = find_50_percent(forward);
    size_t zerophase_50 = find_50_percent(zero_phase);

    std::cout << "Step input at sample: " << step_idx << "\n";
    std::cout << "Forward filter 50% at: " << forward_50
              << " (delay: " << (int)forward_50 - (int)step_idx << ")\n";
    std::cout << "Zero-phase filter 50% at: " << zerophase_50
              << " (delay: " << (int)zerophase_50 - (int)step_idx << ")\n";

    if (std::abs((int)zerophase_50 - (int)step_idx)
        < std::abs((int)forward_50 - (int)step_idx))
    {
        std::cout << "✓ Zero-phase filter has less delay\n";
    }
}

// ============================================================================
// Test 6: Batch Processing (Matrix Filtering)
// ============================================================================

void test6_batch_processing()
{
    std::cout << "\n========================================\n";
    std::cout << "Test 6: Batch Processing\n";
    std::cout << "========================================\n\n";

    const size_t n_samples = 500;
    const size_t n_signals = 5;
    const double fs = 100.0;

    // Create matrix with different frequency signals
    matrixd signals(n_samples, n_signals);
    std::vector<double> signal_freqs = {5, 10, 15, 20, 25};

    for (size_t j = 0; j < n_signals; ++j)
    {
        for (size_t i = 0; i < n_samples; ++i)
        {
            double t = i / fs;
            signals(i, j) =
                std::sin(2.0 * std::numbers::pi * signal_freqs[j] * t);
        }
    }

    // Apply lowpass filter (12 Hz cutoff)
    auto coeffs = signal::butterworth_lowpass(4, 0.24); // 12 / 50
    auto filtered = signal::filtfilt_columns(signals, coeffs);

    std::cout << "Input: " << n_signals << " signals at ";
    for (auto f : signal_freqs)
        std::cout << f << " ";
    std::cout << "Hz\n";
    std::cout << "Filter cutoff: 12 Hz\n\n";

    // Check attenuation
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
        std::cout << "Signal " << j << " (" << signal_freqs[j] << " Hz): "
                  << "attenuation = " << std::setprecision(3) << attenuation
                  << (attenuation > 0.5 ? " (passed)" : " (blocked)") << "\n";
    }
}

int test_filter_matlab()
{
    std::cout << "==============================================\n";
    std::cout << "MSL Validation Program\n";
    std::cout << "Comparing with MATLAB Reference Results\n";
    std::cout << "==============================================\n";

    try
    {
        // Filter Tests
        test1_butterworth_coefficients();
        test2_lowpass_filtering();
        test3_highpass_filtering();
        test4_bandpass_filtering();
        test5_phase_distortion();
        test6_batch_processing();

        std::cout << "\n==============================================\n";
        std::cout << "Validation Complete!\n";
        std::cout << "==============================================\n\n";

        std::cout << "Note: If MATLAB reference files are not found,\n";
        std::cout << "      tests will still run and show MSL results.\n";
        std::cout << "      Run the MATLAB script first to generate\n";
        std::cout << "      reference data for full validation.\n";
    }
    catch (const std::exception &e)
    {
        std::cerr << "\nError: " << e.what() << "\n";
        return 1;
    }

    return 0;
}