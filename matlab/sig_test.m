clc;clear;close all

ori_data = importdata("..\resource\KunmingSSJY.txt");

sig1 = ori_data(:,1);
sig2 = ori_data(:,2);

fs = 50;

pxx = pwelch(sig1);
pyy = pwelch(sig2);
pxy = cpsd(sig1,sig2);