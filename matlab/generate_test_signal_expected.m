clc; clear;

out_dir = fullfile('..', 'test_result', 'matlab', 'signal');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

x = [1.0; 2.0; 3.0; 4.0];
X = fft(x);
restored = ifft(X, 'symmetric');
freqs = (0:3)' .* (8.0 / 4.0);
magnitude = abs(X);

hann_window = hann(5, 'symmetric');
hamming_window = hamming(5, 'symmetric');

b = [0.5, 0.5];
a = 1.0;
filtered = filter(b, a, [1.0; 2.0; 4.0; 8.0]);
filtered_matrix = filter(b, a, [[1.0; 2.0; 4.0; 8.0], ...
                                [2.0; 4.0; 8.0; 16.0]]);

writematrix([real(X), imag(X)], ...
    fullfile(out_dir, 'fft_real_imag.txt'), 'Delimiter', ' ');
writematrix(restored, fullfile(out_dir, 'ifft_real.txt'), 'Delimiter', ' ');
writematrix(freqs, fullfile(out_dir, 'fft_frequencies.txt'), 'Delimiter', ' ');
writematrix(magnitude, ...
    fullfile(out_dir, 'magnitude_spectrum.txt'), 'Delimiter', ' ');
writematrix(hann_window, ...
    fullfile(out_dir, 'hann_window.txt'), 'Delimiter', ' ');
writematrix(hamming_window, ...
    fullfile(out_dir, 'hamming_window.txt'), 'Delimiter', ' ');
writematrix(filtered, fullfile(out_dir, 'filter.txt'), 'Delimiter', ' ');
writematrix(filtered_matrix, ...
    fullfile(out_dir, 'filter_matrix.txt'), 'Delimiter', ' ');
