%% Step 1: 读取 MATLAB 与 C++ 数据
X_mat = readmatrix('matlab_fft.txt');
X_cpp = readmatrix('cpp_fft.txt');

x_ifft_mat = readmatrix('matlab_ifft.txt');
x_ifft_cpp = readmatrix('cpp_ifft.txt');

y_filter_mat = readmatrix('matlab_filter.txt');
y_filter_cpp = readmatrix('cpp_filter.txt');

y_filtfilt_mat = readmatrix('matlab_filtfilt.txt');
y_filtfilt_cpp = readmatrix('cpp_filtfilt.txt');

%% Step 2: 误差计算
X_mat_c = X_mat(:,1) + 1i*X_mat(:,2);
X_cpp_c = X_cpp(:,1) + 1i*X_cpp(:,2);

fprintf('FFT 实部误差: %.3e\n', norm(real(X_mat_c)-real(X_cpp_c))/norm(X_mat_c));
fprintf('FFT 虚部误差: %.3e\n', norm(imag(X_mat_c)-imag(X_cpp_c))/norm(X_mat_c));
fprintf('IFFT 误差: %.3e\n', norm(x_ifft_mat - x_ifft_cpp)/norm(x_ifft_mat));
fprintf('Filter 误差: %.3e\n', norm(y_filter_mat - y_filter_cpp)/norm(y_filter_mat));
fprintf('Filtfilt 误差: %.3e\n', norm(y_filtfilt_mat - y_filtfilt_cpp)/norm(y_filtfilt_mat));

%% Step 3: 可视化比较
figure;
subplot(3,1,1);
plot(real(X_mat_c)); hold on; plot(real(X_cpp_c), '--');
title('FFT 实部比较'); legend('MATLAB','C++');

subplot(3,1,2);
plot(y_filter_mat); hold on; plot(y_filter_cpp, '--');
title('Filter 输出比较'); legend('MATLAB','C++');

subplot(3,1,3);
plot(y_filtfilt_mat); hold on; plot(y_filtfilt_cpp, '--');
title('Filtfilt 输出比较'); legend('MATLAB','C++');
