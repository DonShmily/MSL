%% FourierBandpassFilter - Apply a bandpass filter to a signal using Fourier transform
function osig = FourierBandpassFilter(sig, fs, low, high)
% FourierBandpassFilter - Apply a bandpass filter to a signal using Fourier transform
%
% Syntax: osig = FourierFilter(sig, fs, low, high, order)
%
% Inputs:
%    sig - Input signal (1D array)
%    fs - Sampling frequency (Hz)
%    low - Low cutoff frequency (Hz)
%    high - High cutoff frequency (Hz)
%
% Outputs:
%    osig - Filtered output signal (1D array)
%

sig_fft = fft(sig);
N = length(sig);

% Create a bandpass filter in the frequency domain
H = zeros(size(sig_fft));
H(1:round(low*N/fs)) = 0; % Low cutoff
H(round(low*N/fs)+1:round(high*N/fs)) = 1; % Pass band
H(round(high*N/fs)+1:end) = 0; % High cutoff
% H = H(1:N/2+1); % Keep only positive frequencies
% H = [H, conj(H(end-1:-1:2))]; % Mirror for negative frequencies

osig_fft = sig_fft .* H; % Apply filter in frequency domain
osig = ifft(osig_fft, 'symmetric'); % Inverse FFT to get time domain signal

osig = osig(1:N); % Keep the same length as input signal
end