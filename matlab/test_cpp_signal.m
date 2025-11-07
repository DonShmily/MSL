clc;clear;close all
% 2025.10.23测试通过

ori_data = importdata("..\resource\KunmingSSJY.txt");

fft_filt = FourierBandpassFilter(ori_data(:,1),50,0.25,20);
[b,a] = butter(4,[0.25,20]/(50/2),'bandpass');
butter_filt = filtfilt(b,a,ori_data(:,1));
butter_filt_2 = filtfilt(b,a,ori_data);

fft_filt_cpp = importdata("..\resource\test_result\signal\signal_fft_bandpass.txt");
butter_filt_cpp = importdata("..\resource\test_result\signal\signal_butterworth_bandpass.txt");
butter_filt_2_cpp = importdata("..\resource\test_result\signal\signal_butterworth_bandpass_matrix.txt");

figure
hold on
plot(fft_filt)
plot(fft_filt_cpp)

figure
hold on
plot(butter_filt)
plot(butter_filt_cpp)

%% 
idx = 4;
figure
hold on
plot(butter_filt_2(:,idx))
plot(butter_filt_2_cpp(:,idx))