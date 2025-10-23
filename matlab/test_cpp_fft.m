clc;clear;close all
% 2025.10.23测试通过

ori_data = importdata("..\resource\KunmingSSJY.txt");

fft_res = fft(ori_data(:,1));
fft_fft_res = fft(fft_res);
fft_ifft_res = ifft(fft_res);
fft_fft_ifft_res = ifft(fft_ifft_res);
fft_2d_res = fft(ori_data);

fft_res_cpp = importdata("..\resource\test_result\fft\fft_result.txt");
fft_fft_res_cpp = importdata("..\resource\test_result\fft\fft_result_c2c.txt");
fft_ifft_res_cpp = importdata("..\resource\test_result\fft\ifft_result.txt");
fft_fft_ifft_res_cpp = importdata("..\resource\test_result\fft\ifft_result_c2c.txt");
fft_2d_res_cpp = importdata("..\resource\test_result\fft\fft_matrix_result.txt");
