clc;clear;close all

ori_data = importdata("..\resource\KunmingSSJY.txt");

fft_filt = FourierBandpassFilter(ori_data(:,1),50,0.1,10);
[b,a] = butter(4,[0.1,10]/(50/2),'bandpass');
butter_filt = filtfilt(b,a,ori_data(:,1));

data1 = importdata("..\resource\test_result\signal\signal_butterworth_bandpass.txt");
data2 = importdata("..\\resource\\test_result\\numerical_old\\filtfilt_filtered_data.txt");

figure
hold on
plot(butter_filt)
plot(data1)
plot(data2)