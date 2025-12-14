clc;clear;close all
% 2025.10.23测试通过

ori_data = importdata("..\resource\KunmingSSJY.txt");

len = length(ori_data);
t = (1:len)/50;

co = [1e-2 1e-4 1e-6];
err = co(1) + co(2)*t + co(3)*t.^2;

data1 = ori_data(:,1) + err';

figure
hold on
plot(data1)

data2 = detrend(data1,2);

plot(data2)

p = polyfit(t, data1, 2);
data3 = data1 - polyval(p,t)';
plot(data3)

data4 = importdata("..\resource\test_result\signal\detrended_data.txt");
plot(data4)