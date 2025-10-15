clc;clear;close all

%% Step 1: 生成测试信号
% fs = 1000;                % 采样率
% t = (0:1/fs:1-1/fs)';     % 时间轴
% x = sin(2*pi*50*t) + 0.5*sin(2*pi*120*t) + 0.8*cos(2*pi*80*t); % 混合信号
% 
% writematrix([t x], 'test_signal.txt', 'Delimiter', ' ');
x = importdata("test_signal.txt");
x = x(:,2);

%% Step 2: MATLAB FFT / IFFT
X = fft(x);
x_ifft = ifft(X, 'symmetric');

% 分开实部与虚部保存
writematrix([real(X) imag(X)], 'matlab_fft.txt', 'Delimiter', ' ');
writematrix(x_ifft, 'matlab_ifft.txt', 'Delimiter', ' ');

%% Step 3: 设计滤波器并滤波
order = 4;
fc = 0.5; % 归一化截止频率（0~1, 1 对应 Nyquist）
[b, a] = butter(4, [0.05 0.5], 'bandpass');
% [b, a] = butter(order, 0.02, 'high');
disp(order)
disp(b)
disp(a)

y_filter = filter(b, a, x);
y_filtfilt = filtfilt(b, a, x);

% 保存滤波器系数与输出
writematrix([b(:) a(:)], 'filter_coeffs.txt', 'Delimiter', ' ');
writematrix(y_filter, 'matlab_filter.txt', 'Delimiter', ' ');
writematrix(y_filtfilt, 'matlab_filtfilt.txt', 'Delimiter', ' ');
