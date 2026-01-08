clc;clear;

ori_data = importdata("../resource/KunmingSSJY.txt");
acc1 = ori_data(:,1);

fs = 50;
low = 0.1;
high = 10;

[b,a] = butter(4,[low high]/(fs/2));
filtfilt_data = filtfilt(b,a,acc1);

cpp_data = importdata("..\\resource\\test_result\\numerical_old\\filtfilt_filtered_data.txt");

figure;hold on;
% plot(acc1,'k');
plot(filtfilt_data);
plot(cpp_data);