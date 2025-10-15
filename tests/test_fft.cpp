/**
 * Examples and tests for MSL FFT module
 */

#include "test.hpp"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <vector>

#include "fft/fft.hpp"

#include "helper.hpp"

using namespace msl;

// ============================================================================
// Test 1: Basic 1D FFT of sine wave
// ============================================================================

void test_1d_fft_sine()
{
    std::cout << "========================================\n";
    std::cout << "Test 1: 1D FFT of Sine Wave\n";
    std::cout << "========================================\n\n";

    // Generate sine wave: f = 5 Hz, sampling rate = 100 Hz
    const size_t n = 128;
    const double fs = 100.0; // Sampling rate
    const double f0 = 5.0;   // Signal frequency

    std::vector<double> signal(n);
    for (size_t i = 0; i < n; ++i)
    {
        double t = i / fs;
        signal[i] = std::sin(2.0 * std::numbers::pi * f0 * t);
    }

    print_vector(signal, "Input signal (first 10)", 10);

    // Compute FFT
    auto fft_result = fft::fft(signal);
    print_complex_vector(fft_result, "FFT result (first 10)", 10);

    // Compute magnitude spectrum
    auto magnitude = fft::magnitude_spectrum(fft_result);
    print_vector(magnitude, "Magnitude spectrum (first 10)", 10);

    // Find peak frequency
    size_t peak_idx = 0;
    double peak_val = 0.0;
    for (size_t i = 1; i < magnitude.size(); ++i)
    { // Skip DC
        if (magnitude[i] > peak_val)
        {
            peak_val = magnitude[i];
            peak_idx = i;
        }
    }

    auto freqs = fft::fft_frequencies(n, fs);
    std::cout << "\nPeak at index " << peak_idx
              << ", frequency = " << freqs[peak_idx] << " Hz\n";
    std::cout << "Expected frequency: " << f0 << " Hz\n\n";

    // Test inverse FFT
    auto reconstructed = fft::ifft_real(fft_result, n);
    print_vector(reconstructed, "Reconstructed signal (first 10)", 10);

    // Compute reconstruction error
    double max_error = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        max_error = std::max(max_error, std::abs(signal[i] - reconstructed[i]));
    }
    std::cout << "Reconstruction error: " << max_error << "\n\n";
}

// ============================================================================
// Test 2: Multi-frequency signal
// ============================================================================

void test_multi_frequency()
{
    std::cout << "========================================\n";
    std::cout << "Test 2: Multi-Frequency Signal\n";
    std::cout << "========================================\n\n";

    // Signal: 3 Hz + 7 Hz + 15 Hz
    const size_t n = 256;
    const double fs = 100.0;

    std::vector<double> signal(n);
    for (size_t i = 0; i < n; ++i)
    {
        double t = i / fs;
        signal[i] = std::sin(2.0 * std::numbers::pi * 3.0 * t)
                    + 0.5 * std::sin(2.0 * std::numbers::pi * 7.0 * t)
                    + 0.3 * std::sin(2.0 * std::numbers::pi * 15.0 * t);
    }

    // Apply Hamming window
    std::vector<double> windowed = signal;
    fft::apply_hamming_window(windowed);

    // Compute FFT
    auto fft_result = fft::fft(windowed);
    auto magnitude = fft::magnitude_spectrum(fft_result);
    auto freqs = fft::fft_frequencies(n, fs);

    // Find peaks
    std::cout << "Detected frequency peaks:\n";
    for (size_t i = 1; i < magnitude.size() - 1; ++i)
    {
        if (magnitude[i] > magnitude[i - 1] && magnitude[i] > magnitude[i + 1]
            && magnitude[i] > 10.0)
        { // Threshold
            std::cout << "  Frequency: " << freqs[i] << " Hz, "
                      << "Magnitude: " << magnitude[i] << "\n";
        }
    }
    std::cout << "\nExpected frequencies: 3, 7, 15 Hz\n\n";
}

// ============================================================================
// Test 3: 2D FFT on matrix (column-wise)
// ============================================================================

void test_2d_fft_columns()
{
    std::cout << "========================================\n";
    std::cout << "Test 3: 2D FFT (Column-wise)\n";
    std::cout << "========================================\n\n";

    // Create matrix with different sine waves in each column
    const size_t n_rows = 64;
    const size_t n_cols = 3;
    const double fs = 100.0;

    matrixd signals(n_rows, n_cols);

    // Column 0: 5 Hz
    // Column 1: 10 Hz
    // Column 2: 15 Hz
    std::vector<double> frequencies = {5.0, 10.0, 15.0};

    for (size_t j = 0; j < n_cols; ++j)
    {
        for (size_t i = 0; i < n_rows; ++i)
        {
            double t = i / fs;
            signals(i, j) =
                std::sin(2.0 * std::numbers::pi * frequencies[j] * t);
        }
    }

    std::cout << "Input matrix: " << n_rows << "x" << n_cols << "\n";
    std::cout << "Each column contains a sine wave at different frequency\n\n";

    // Compute FFT on each column
    auto fft_result = fft::fft_columns(signals);

    std::cout << "FFT result matrix: " << fft_result.rows() << "x"
              << fft_result.cols() << "\n\n";

    // Find peak in each column
    auto freqs = fft::fft_frequencies(n_rows, fs);

    for (size_t j = 0; j < n_cols; ++j)
    {
        size_t peak_idx = 0;
        double peak_mag = 0.0;

        for (size_t i = 1; i < fft_result.rows(); ++i)
        {
            double mag = std::abs(fft_result(i, j));
            if (mag > peak_mag)
            {
                peak_mag = mag;
                peak_idx = i;
            }
        }

        std::cout << "Column " << j << " peak: " << freqs[peak_idx] << " Hz "
                  << "(expected: " << frequencies[j] << " Hz)\n";
    }
    std::cout << "\n";

    // Test inverse FFT
    auto reconstructed = fft::ifft_columns_real(fft_result, n_rows);

    std::cout << "Reconstructed matrix: " << reconstructed.rows() << "x"
              << reconstructed.cols() << "\n";

    // Check error
    double max_error = 0.0;
    for (size_t i = 0; i < n_rows; ++i)
    {
        for (size_t j = 0; j < n_cols; ++j)
        {
            max_error = std::max(max_error,
                                 std::abs(signals(i, j) - reconstructed(i, j)));
        }
    }
    std::cout << "Reconstruction error: " << max_error << "\n\n";
}

// ============================================================================
// Test 4: Window functions comparison
// ============================================================================

void test_window_functions()
{
    std::cout << "========================================\n";
    std::cout << "Test 5: Window Functions\n";
    std::cout << "========================================\n\n";

    const size_t n = 64;
    const double fs = 100.0;
    const double f0 = 10.0;

    // Generate signal
    std::vector<double> signal(n);
    for (size_t i = 0; i < n; ++i)
    {
        double t = i / fs;
        signal[i] = std::sin(2.0 * std::numbers::pi * f0 * t);
    }

    // Test without window
    {
        auto fft_result = fft::fft(signal);
        auto magnitude = fft::magnitude_spectrum(fft_result);

        double spectral_leakage = 0.0;
        for (size_t i = 0; i < magnitude.size(); ++i)
        {
            if (i != 6 && i != 7)
            { // Skip main peak bins
                spectral_leakage += magnitude[i];
            }
        }

        std::cout << "No window - Spectral leakage: " << spectral_leakage
                  << "\n";
    }

    // Test with Hamming window
    {
        std::vector<double> windowed = signal;
        fft::apply_hamming_window(windowed);

        auto fft_result = fft::fft(windowed);
        auto magnitude = fft::magnitude_spectrum(fft_result);

        double spectral_leakage = 0.0;
        for (size_t i = 0; i < magnitude.size(); ++i)
        {
            if (i != 6 && i != 7)
            {
                spectral_leakage += magnitude[i];
            }
        }

        std::cout << "Hamming window - Spectral leakage: " << spectral_leakage
                  << "\n";
    }

    // Test with Hann window
    {
        std::vector<double> windowed = signal;
        fft::apply_hann_window(windowed);

        auto fft_result = fft::fft(windowed);
        auto magnitude = fft::magnitude_spectrum(fft_result);

        double spectral_leakage = 0.0;
        for (size_t i = 0; i < magnitude.size(); ++i)
        {
            if (i != 6 && i != 7)
            {
                spectral_leakage += magnitude[i];
            }
        }

        std::cout << "Hann window - Spectral leakage: " << spectral_leakage
                  << "\n";
    }

    std::cout << "\n";
}

// ============================================================================
// Test 5: Power spectrum analysis
// ============================================================================

void test_power_spectrum()
{
    std::cout << "========================================\n";
    std::cout << "Test 6: Power Spectrum Analysis\n";
    std::cout << "========================================\n\n";

    const size_t n = 256;
    const double fs = 1000.0; // 1 kHz sampling

    // Generate signal with noise
    std::vector<double> signal(n);
    for (size_t i = 0; i < n; ++i)
    {
        double t = i / fs;
        // Signal: 50 Hz and 120 Hz components + noise
        signal[i] = std::sin(2.0 * std::numbers::pi * 50.0 * t)
                    + 0.5 * std::sin(2.0 * std::numbers::pi * 120.0 * t);
        // Add small noise (simplified - normally use random)
        signal[i] += 0.1 * std::sin(2.0 * std::numbers::pi * 327.0 * t);
    }

    // Apply window and compute FFT
    fft::apply_hamming_window(signal);
    auto fft_result = fft::fft(signal);

    // Compute power spectrum
    auto power = fft::power_spectrum(fft_result);
    auto freqs = fft::fft_frequencies(n, fs);

    // Find peaks in power spectrum
    std::cout << "Significant frequency components (power > 1000):\n";
    for (size_t i = 1; i < power.size() - 1; ++i)
    {
        if (power[i] > 1000.0 && power[i] > power[i - 1]
            && power[i] > power[i + 1])
        {
            std::cout << "  f = " << freqs[i] << " Hz, "
                      << "Power = " << power[i] << "\n";
        }
    }

    // Compute total power
    double total_power = 0.0;
    for (auto p : power)
    {
        total_power += p;
    }

    std::cout << "\nTotal power: " << total_power << "\n\n";
}

// ============================================================================
// Test 6: Performance benchmark
// ============================================================================

void test_performance()
{
    std::cout << "========================================\n";
    std::cout << "Test 7: Performance Benchmark\n";
    std::cout << "========================================\n\n";

    std::vector<size_t> sizes = {64, 128, 256, 512, 1024, 2048, 4096};

    std::cout << std::setw(10) << "Size" << std::setw(15) << "Time (ms)"
              << std::setw(20) << "Time per sample (μs)\n";
    std::cout << std::string(45, '-') << "\n";

    for (auto n : sizes)
    {
        // Generate signal
        std::vector<double> signal(n);
        for (size_t i = 0; i < n; ++i)
        {
            signal[i] = std::sin(2.0 * std::numbers::pi * i / n);
        }

        // Benchmark
        const int iterations = 1000;
        auto start = std::chrono::high_resolution_clock::now();

        for (int iter = 0; iter < iterations; ++iter)
        {
            auto result = fft::fft(signal);
            // Prevent optimization
            volatile double dummy = result[0].real();
            (void)dummy;
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration =
            std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        double avg_ms = duration.count() / 1000.0 / iterations;
        double per_sample_us =
            duration.count() / static_cast<double>(iterations) / n;

        std::cout << std::setw(10) << n << std::setw(15) << std::fixed
                  << std::setprecision(3) << avg_ms << std::setw(20)
                  << std::fixed << std::setprecision(6) << per_sample_us
                  << "\n";
    }
    std::cout << "\n";
}

// ============================================================================
// Test 7: Complex input FFT
// ============================================================================

void test_complex_fft()
{
    std::cout << "========================================\n";
    std::cout << "Test 8: Complex Input FFT\n";
    std::cout << "========================================\n\n";

    const size_t n = 64;

    // Generate complex signal
    std::vector<std::complex<double>> signal(n);
    for (size_t i = 0; i < n; ++i)
    {
        double t = i / 64.0;
        // Complex exponential: e^(i*2π*5*t)
        signal[i] = std::exp(
            std::complex<double>(0.0, 2.0 * std::numbers::pi * 5.0 * t));
    }

    print_complex_vector(signal, "Complex input signal (first 10)", 10);

    // Compute FFT
    auto fft_result = fft::fft(signal);

    std::cout << "\nFFT magnitude at each frequency:\n";
    for (size_t i = 0; i < std::min(size_t(15), fft_result.size()); ++i)
    {
        std::cout << "  f[" << i << "] = " << std::abs(fft_result[i]) << "\n";
    }

    std::cout << "\nExpected: peak at frequency bin 5\n\n";

    // Inverse FFT
    auto reconstructed = fft::ifft(fft_result);

    // Check error
    double max_error = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        max_error = std::max(max_error, std::abs(signal[i] - reconstructed[i]));
    }
    std::cout << "Reconstruction error: " << max_error << "\n\n";
}

// ============================================================================
// Main
// ============================================================================

int test_fft()
{
    try
    {
        test_1d_fft_sine();
        test_multi_frequency();
        test_2d_fft_columns();
        test_window_functions();
        test_power_spectrum();
        test_complex_fft();
        test_performance();

        std::cout << "========================================\n";
        std::cout << "All FFT tests completed successfully!\n";
        std::cout << "========================================\n";
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}