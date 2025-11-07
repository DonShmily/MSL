clc;clear;close all
% 2025.11.08测试通过

ori_data = importdata("..\resource\KunmingSSJY.txt");

pxx = pwelch(ori_data(:,1),hann(1024),512,1024);
pxy = cpsd(ori_data(:,1),ori_data(:,2),hann(1024),512,1024);
pxy_mag = abs(pxy);

psd_cpp = importdata("..\resource\test_result\signal\psd_1d.txt");
cpsd_cpp = importdata("..\resource\test_result\signal\cpsd_2d.txt");
cpsd_cpp_mag =importdata("..\resource\test_result\signal\cpsd_2d_mag.txt");

figure
hold on
plot(pxx/max(pxx))
plot(psd_cpp/max(psd_cpp))

figure
hold on
plot(pxy_mag/max(pxy_mag))
plot(cpsd_cpp_mag/max(cpsd_cpp_mag))
