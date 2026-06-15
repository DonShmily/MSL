script_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(script_dir, '..', '..', '..');
data_dir = fullfile(repo_root, 'test_result', 'signal', 'matlab_compare');

fft_file = fullfile(data_dir, 'signal_fft.txt');
window_file = fullfile(data_dir, 'signal_windows.txt');
filter_file = fullfile(data_dir, 'signal_filter.txt');
if ~isfile(fft_file) || ~isfile(window_file) || ~isfile(filter_file)
    error('Missing comparison data. Run "xmake run test_signal" first.');
end

fft_data = readmatrix(fft_file);
x = fft_data(:, 1);
X_msl = fft_data(:, 2) + 1i * fft_data(:, 3);
restored = fft_data(:, 4);
X_ref = fft(x);

assert(max(abs(X_msl - X_ref)) < 1e-12, ...
    'FFT comparison failed.');
assert(max(abs(restored - real(ifft(X_ref)))) < 1e-12, ...
    'IFFT comparison failed.');

windows = readmatrix(window_file);
n = size(windows, 1);
k = (0:n-1)';
hann_ref = 0.5 - 0.5 * cos(2 * pi * k / (n - 1));
hamming_ref = 0.54 - 0.46 * cos(2 * pi * k / (n - 1));

assert(max(abs(windows(:, 1) - hann_ref)) < 1e-12, ...
    'Hann window comparison failed.');
assert(max(abs(windows(:, 2) - hamming_ref)) < 1e-12, ...
    'Hamming window comparison failed.');

filter_data = readmatrix(filter_file);
filter_ref = filter([0.5, 0.5], 1.0, filter_data(:, 1));
assert(max(abs(filter_data(:, 2) - filter_ref)) < 1e-12, ...
    'Filter comparison failed.');

fprintf('Signal validation passed. FFT error = %.3e\n', ...
    max(abs(X_msl - X_ref)));
