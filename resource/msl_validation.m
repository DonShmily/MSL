%% MSL Library Validation Script
% This script generates test data and expected results for validating
% MSL FFT and Filter implementations against MATLAB

clear; clc; close all;

%% Configuration
output_dir = 'msl_validation_data';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('==============================================\n');
fprintf('MSL Library Validation Test Data Generation\n');
fprintf('==============================================\n\n');

%% ====================================================================
%% Test 1: FFT - Single Frequency Sine Wave
%% ====================================================================
fprintf('Test 1: FFT - Single Frequency Sine Wave\n');
fprintf('------------------------------------------\n');

% Parameters
n1 = 128;
fs1 = 100.0;  % Sampling rate
f0_1 = 5.0;   % Signal frequency

% Generate signal
t1 = (0:n1-1) / fs1;
signal1 = sin(2*pi*f0_1*t1)';

% Compute FFT
fft_result1 = fft(signal1);
fft_result1_half = fft_result1(1:floor(n1/2)+1);  % Only positive frequencies

% Frequency axis
freqs1 = (0:floor(n1/2)) * fs1 / n1;

% Magnitude and phase
magnitude1 = abs(fft_result1_half);
phase1 = angle(fft_result1_half);

% Save results
save(fullfile(output_dir, 'test1_fft_sine.mat'), ...
    'signal1', 'fft_result1_half', 'magnitude1', 'phase1', 'freqs1', ...
    'n1', 'fs1', 'f0_1');

fprintf('  Signal length: %d\n', n1);
fprintf('  Sampling rate: %.1f Hz\n', fs1);
fprintf('  Signal frequency: %.1f Hz\n', f0_1);
fprintf('  Peak magnitude: %.4f at %.2f Hz\n', max(magnitude1), freqs1(magnitude1 == max(magnitude1)));
fprintf('  Saved to: test1_fft_sine.mat\n\n');

%% ====================================================================
%% Test 2: FFT - Multi-frequency Signal
%% ====================================================================
fprintf('Test 2: FFT - Multi-frequency Signal\n');
fprintf('--------------------------------------\n');

% Parameters
n2 = 256;
fs2 = 100.0;
freqs_signal2 = [3.0, 7.0, 15.0];  % Three frequencies
amps2 = [1.0, 0.5, 0.3];

% Generate signal
t2 = (0:n2-1) / fs2;
signal2 = zeros(n2, 1);
for i = 1:length(freqs_signal2)
    signal2 = signal2 + amps2(i) * sin(2*pi*freqs_signal2(i)*t2)';
end

% Apply Hamming window
window2 = hamming(n2);
signal2_windowed = signal2 .* window2;

% Compute FFT
fft_result2 = fft(signal2_windowed);
fft_result2_half = fft_result2(1:floor(n2/2)+1);

% Frequency axis
freqs2 = (0:floor(n2/2)) * fs2 / n2;

% Magnitude
magnitude2 = abs(fft_result2_half);

% Save results
save(fullfile(output_dir, 'test2_fft_multifreq.mat'), ...
    'signal2', 'signal2_windowed', 'fft_result2_half', 'magnitude2', 'freqs2', ...
    'n2', 'fs2', 'freqs_signal2', 'amps2');

fprintf('  Signal length: %d\n', n2);
fprintf('  Frequencies: %.1f, %.1f, %.1f Hz\n', freqs_signal2);
fprintf('  Window: Hamming\n');
fprintf('  Saved to: test2_fft_multifreq.mat\n\n');

%% ====================================================================
%% Test 3: FFT - 2D Matrix (Column-wise)
%% ====================================================================
fprintf('Test 3: FFT - 2D Matrix (Column-wise)\n');
fprintf('---------------------------------------\n');

% Parameters
n_rows3 = 64;
n_cols3 = 3;
fs3 = 100.0;
col_freqs3 = [5.0, 10.0, 15.0];  % Each column has different frequency

% Generate matrix
t3 = (0:n_rows3-1) / fs3;
signals3 = zeros(n_rows3, n_cols3);
for j = 1:n_cols3
    signals3(:, j) = sin(2*pi*col_freqs3(j)*t3)';
end

% Compute FFT for each column
fft_matrix3 = zeros(floor(n_rows3/2)+1, n_cols3);
for j = 1:n_cols3
    fft_col = fft(signals3(:, j));
    fft_matrix3(:, j) = fft_col(1:floor(n_rows3/2)+1);
end

% Magnitude
magnitude_matrix3 = abs(fft_matrix3);

% Frequency axis
freqs3 = (0:floor(n_rows3/2)) * fs3 / n_rows3;

% Save results
save(fullfile(output_dir, 'test3_fft_matrix.mat'), ...
    'signals3', 'fft_matrix3', 'magnitude_matrix3', 'freqs3', ...
    'n_rows3', 'n_cols3', 'fs3', 'col_freqs3');

fprintf('  Matrix size: %dx%d\n', n_rows3, n_cols3);
fprintf('  Column frequencies: %.1f, %.1f, %.1f Hz\n', col_freqs3);
fprintf('  Saved to: test3_fft_matrix.mat\n\n');

%% ====================================================================
%% Test 4: IFFT - Reconstruction
%% ====================================================================
fprintf('Test 4: IFFT - Reconstruction\n');
fprintf('------------------------------\n');

% Use signal from Test 1
fft_full4 = fft(signal1);
reconstructed4 = ifft(fft_full4);
reconstructed4 = real(reconstructed4);  % Take real part

% Error
error4 = max(abs(signal1 - reconstructed4));

% Save results
save(fullfile(output_dir, 'test4_ifft_reconstruction.mat'), ...
    'signal1', 'fft_full4', 'reconstructed4', 'error4');

fprintf('  Reconstruction error: %.2e\n', error4);
fprintf('  Saved to: test4_ifft_reconstruction.mat\n\n');

%% ====================================================================
%% Test 5: Butterworth Lowpass Filter Design
%% ====================================================================
fprintf('Test 5: Butterworth Lowpass Filter\n');
fprintf('-----------------------------------\n');

% Parameters
order5 = 4;
fc5_hz = 30.0;    % Cutoff frequency in Hz
fs5 = 1000.0;     % Sampling rate
fc5_norm = fc5_hz / (fs5/2);  % Normalized frequency

% Design filter
[b5, a5] = butter(order5, fc5_norm, 'low');

% Generate test signal (10 Hz + 100 Hz noise)
n5 = 1000;
t5 = (0:n5-1) / fs5;
signal5 = sin(2*pi*10*t5)' + 0.5*sin(2*pi*100*t5)';

% Apply filter (single-pass)
filtered5_single = filter(b5, a5, signal5);

% Apply filter (zero-phase)
filtered5_zerophase = filtfilt(b5, a5, signal5);

% Frequency response
[H5, f5] = freqz(b5, a5, 512, fs5);
magnitude_response5 = abs(H5);
phase_response5 = angle(H5);

% Save results
save(fullfile(output_dir, 'test5_butter_lowpass.mat'), ...
    'b5', 'a5', 'signal5', 'filtered5_single', 'filtered5_zerophase', ...
    'H5', 'f5', 'magnitude_response5', 'phase_response5', ...
    'order5', 'fc5_hz', 'fs5', 'fc5_norm');

fprintf('  Filter order: %d\n', order5);
fprintf('  Cutoff frequency: %.1f Hz (normalized: %.4f)\n', fc5_hz, fc5_norm);
fprintf('  Sampling rate: %.1f Hz\n', fs5);
fprintf('  Saved to: test5_butter_lowpass.mat\n\n');

%% ====================================================================
%% Test 6: Butterworth Highpass Filter Design
%% ====================================================================
fprintf('Test 6: Butterworth Highpass Filter\n');
fprintf('------------------------------------\n');

% Parameters
order6 = 4;
fc6_hz = 1.0;
fs6 = 100.0;
fc6_norm = fc6_hz / (fs6/2);

% Design filter
[b6, a6] = butter(order6, fc6_norm, 'high');

% Generate test signal with DC offset
n6 = 500;
t6 = (0:n6-1) / fs6;
dc_offset6 = 5.0;
signal6 = dc_offset6 + sin(2*pi*2*t6)';

% Apply filter
filtered6_zerophase = filtfilt(b6, a6, signal6);

% Frequency response
[H6, f6] = freqz(b6, a6, 512, fs6);

% Save results
save(fullfile(output_dir, 'test6_butter_highpass.mat'), ...
    'b6', 'a6', 'signal6', 'filtered6_zerophase', 'H6', 'f6', ...
    'order6', 'fc6_hz', 'fs6', 'fc6_norm', 'dc_offset6');

fprintf('  Filter order: %d\n', order6);
fprintf('  Cutoff frequency: %.1f Hz (normalized: %.4f)\n', fc6_hz, fc6_norm);
fprintf('  DC offset: %.1f\n', dc_offset6);
fprintf('  Saved to: test6_butter_highpass.mat\n\n');

%% ====================================================================
%% Test 7: Butterworth Bandpass Filter Design
%% ====================================================================
fprintf('Test 7: Butterworth Bandpass Filter\n');
fprintf('------------------------------------\n');

% Parameters
order7 = 4;
fc_low7_hz = 30.0;
fc_high7_hz = 100.0;
fs7 = 1000.0;
fc_low7_norm = fc_low7_hz / (fs7/2);
fc_high7_norm = fc_high7_hz / (fs7/2);

% Design filter
[b7, a7] = butter(order7, [fc_low7_norm, fc_high7_norm], 'bandpass');

% Generate test signal (10, 50, 200 Hz)
n7 = 1024;
t7 = (0:n7-1) / fs7;
signal7 = sin(2*pi*10*t7)' + sin(2*pi*50*t7)' + sin(2*pi*200*t7)';

% Apply filter
filtered7_zerophase = filtfilt(b7, a7, signal7);

% Frequency response
[H7, f7] = freqz(b7, a7, 512, fs7);

% Save results
save(fullfile(output_dir, 'test7_butter_bandpass.mat'), ...
    'b7', 'a7', 'signal7', 'filtered7_zerophase', 'H7', 'f7', ...
    'order7', 'fc_low7_hz', 'fc_high7_hz', 'fs7', ...
    'fc_low7_norm', 'fc_high7_norm');

fprintf('  Filter order: %d\n', order7);
fprintf('  Passband: %.1f - %.1f Hz\n', fc_low7_hz, fc_high7_hz);
fprintf('  Normalized: %.4f - %.4f\n', fc_low7_norm, fc_high7_norm);
fprintf('  Saved to: test7_butter_bandpass.mat\n\n');

%% ====================================================================
%% Test 8: Filter Coefficients Comparison
%% ====================================================================
fprintf('Test 8: Filter Coefficients Export\n');
fprintf('-----------------------------------\n');

% Various filter configurations for coefficient comparison
configs = struct();

% Config 1: 2nd order lowpass
[b_c1, a_c1] = butter(2, 0.2, 'low');
configs.lowpass_o2_fc02.b = b_c1;
configs.lowpass_o2_fc02.a = a_c1;

% Config 2: 4th order lowpass
[b_c2, a_c2] = butter(4, 0.2, 'low');
configs.lowpass_o4_fc02.b = b_c2;
configs.lowpass_o4_fc02.a = a_c2;

% Config 3: 4th order highpass
[b_c3, a_c3] = butter(4, 0.3, 'high');
configs.highpass_o4_fc03.b = b_c3;
configs.highpass_o4_fc03.a = a_c3;

% Config 4: 2nd order bandpass
[b_c4, a_c4] = butter(2, [0.1, 0.4], 'bandpass');
configs.bandpass_o2_fc01_04.b = b_c4;
configs.bandpass_o2_fc01_04.a = a_c4;

% Config 5: 6th order lowpass
[b_c5, a_c5] = butter(6, 0.15, 'low');
configs.lowpass_o6_fc015.b = b_c5;
configs.lowpass_o6_fc015.a = a_c5;

% Save all configurations
save(fullfile(output_dir, 'test8_filter_coefficients.mat'), 'configs');

fprintf('  Exported %d filter configurations\n', length(fieldnames(configs)));
fprintf('  Saved to: test8_filter_coefficients.mat\n\n');

%% ====================================================================
%% Test 9: Matrix Filtering (Batch Processing)
%% ====================================================================
fprintf('Test 9: Matrix Filtering (Batch)\n');
fprintf('---------------------------------\n');

% Parameters
n_samples9 = 500;
n_signals9 = 5;
fs9 = 100.0;

% Create matrix with signals at different frequencies
signals9 = zeros(n_samples9, n_signals9);
signal_freqs9 = [5, 10, 15, 20, 25];  % Hz

t9 = (0:n_samples9-1)' / fs9;
for j = 1:n_signals9
    signals9(:, j) = sin(2*pi*signal_freqs9(j)*t9);
end

% Design lowpass filter at 12 Hz
fc9_hz = 12.0;
fc9_norm = fc9_hz / (fs9/2);
[b9, a9] = butter(4, fc9_norm, 'low');

% Filter each column
filtered9 = zeros(size(signals9));
for j = 1:n_signals9
    filtered9(:, j) = filtfilt(b9, a9, signals9(:, j));
end

% Save results
save(fullfile(output_dir, 'test9_matrix_filtering.mat'), ...
    'signals9', 'filtered9', 'b9', 'a9', ...
    'n_samples9', 'n_signals9', 'fs9', 'signal_freqs9', 'fc9_hz');

fprintf('  Matrix size: %dx%d\n', n_samples9, n_signals9);
fprintf('  Signal frequencies: '); fprintf('%.0f ', signal_freqs9); fprintf('Hz\n');
fprintf('  Filter cutoff: %.1f Hz\n', fc9_hz);
fprintf('  Saved to: test9_matrix_filtering.mat\n\n');

%% ====================================================================
%% Test 10: Edge Cases
%% ====================================================================
fprintf('Test 10: Edge Cases\n');
fprintf('-------------------\n');

% Small signal
small_signal = [1, 2, 3, 4]';
small_fft = fft(small_signal);

% Single sample (degenerate case)
single_sample = 1.0;
single_fft = fft(single_sample);

% Power of 2 sizes
sizes_pot = [16, 32, 64, 128, 256, 512, 1024];
pot_signals = cell(length(sizes_pot), 1);
pot_ffts = cell(length(sizes_pot), 1);

for i = 1:length(sizes_pot)
    n = sizes_pot(i);
    t = (0:n-1) / 100.0;
    pot_signals{i} = sin(2*pi*5*t)';
    pot_ffts{i} = fft(pot_signals{i});
end

% Save results
save(fullfile(output_dir, 'test10_edge_cases.mat'), ...
    'small_signal', 'small_fft', 'single_sample', 'single_fft', ...
    'sizes_pot', 'pot_signals', 'pot_ffts');

fprintf('  Small signal: length %d\n', length(small_signal));
fprintf('  Power-of-2 sizes: '); fprintf('%d ', sizes_pot); fprintf('\n');
fprintf('  Saved to: test10_edge_cases.mat\n\n');

%% ====================================================================
%% Generate Comparison CSV Files
%% ====================================================================
fprintf('Generating CSV files for easy comparison...\n');
fprintf('-------------------------------------------\n');

% Test 1: FFT magnitudes
csvwrite(fullfile(output_dir, 'test1_fft_magnitude.csv'), ...
    [freqs1', magnitude1]);
fprintf('  Written: test1_fft_magnitude.csv\n');

% Test 5: Filter coefficients
csvwrite(fullfile(output_dir, 'test5_filter_b_coeffs.csv'), b5);
csvwrite(fullfile(output_dir, 'test5_filter_a_coeffs.csv'), a5);
fprintf('  Written: test5_filter_*_coeffs.csv\n');

% Test 5: Filtered signals
csvwrite(fullfile(output_dir, 'test5_original_signal.csv'), signal5);
csvwrite(fullfile(output_dir, 'test5_filtered_zerophase.csv'), filtered5_zerophase);
fprintf('  Written: test5_*_signal.csv\n');

fprintf('\n');

%% ====================================================================
%% Generate Summary Report
%% ====================================================================
fprintf('==============================================\n');
fprintf('Summary Report\n');
fprintf('==============================================\n\n');

summary_file = fullfile(output_dir, 'validation_summary.txt');
fid = fopen(summary_file, 'w');

fprintf(fid, 'MSL Library Validation Test Data\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));

fprintf(fid, 'Test Files Generated:\n');
fprintf(fid, '---------------------\n');
fprintf(fid, '1. test1_fft_sine.mat - Single frequency FFT\n');
fprintf(fid, '2. test2_fft_multifreq.mat - Multi-frequency FFT\n');
fprintf(fid, '3. test3_fft_matrix.mat - Matrix FFT (column-wise)\n');
fprintf(fid, '4. test4_ifft_reconstruction.mat - IFFT reconstruction\n');
fprintf(fid, '5. test5_butter_lowpass.mat - Butterworth lowpass filter\n');
fprintf(fid, '6. test6_butter_highpass.mat - Butterworth highpass filter\n');
fprintf(fid, '7. test7_butter_bandpass.mat - Butterworth bandpass filter\n');
fprintf(fid, '8. test8_filter_coefficients.mat - Various filter configs\n');
fprintf(fid, '9. test9_matrix_filtering.mat - Batch filtering\n');
fprintf(fid, '10. test10_edge_cases.mat - Edge cases\n\n');

fprintf(fid, 'CSV Files for Quick Comparison:\n');
fprintf(fid, '-------------------------------\n');
fprintf(fid, '- test1_fft_magnitude.csv\n');
fprintf(fid, '- test5_filter_b_coeffs.csv\n');
fprintf(fid, '- test5_filter_a_coeffs.csv\n');
fprintf(fid, '- test5_original_signal.csv\n');
fprintf(fid, '- test5_filtered_zerophase.csv\n\n');

fprintf(fid, 'Validation Instructions:\n');
fprintf(fid, '------------------------\n');
fprintf(fid, '1. Run your MSL FFT/Filter tests\n');
fprintf(fid, '2. Compare MSL outputs with corresponding .mat files\n');
fprintf(fid, '3. Check maximum absolute error\n');
fprintf(fid, '4. Acceptable error: < 1e-10 for FFT, < 1e-6 for filters\n\n');

fprintf(fid, 'Key Parameters:\n');
fprintf(fid, '---------------\n');
fprintf(fid, 'FFT Test 1: n=%d, fs=%.1f Hz, f0=%.1f Hz\n', n1, fs1, f0_1);
fprintf(fid, 'Filter Test 5: order=%d, fc=%.1f Hz, fs=%.1f Hz\n', order5, fc5_hz, fs5);
fprintf(fid, 'Filter Test 7: order=%d, passband=%.1f-%.1f Hz\n', order7, fc_low7_hz, fc_high7_hz);

fclose(fid);

fprintf('Generated: validation_summary.txt\n\n');

fprintf('==============================================\n');
fprintf('All test data generated successfully!\n');
fprintf('Output directory: %s\n', output_dir);
fprintf('==============================================\n');