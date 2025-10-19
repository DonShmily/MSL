clc;clear;

ori_data = importdata("sig.txt");

filt_data = FourierBandpassFilter(ori_data, 50, 2.5, 7.5);

data_cpp1 = importdata("fft_filter.txt");
data_cpp2 = importdata("fft_filter_manual_freq_response.txt");

figure
hold on
% plot(ori_data)
plot(filt_data)
plot(data_cpp1)
% plot(data_cpp2)