clc;clear;close all
% 2025.10.23测试通过

ori_data = importdata("..\resource\KunmingSSJY.txt");

intgral_1d = cumtrapz(ori_data(:,1));
intgral_1d_nonuni = cumtrapz(0.1,ori_data(:,1));
intgral_2d = cumtrapz(ori_data);
intgral_2d_nonuni = cumtrapz(0.1,ori_data);

intgral_1d_cpp = importdata("..\resource\test_result\integral\integral_1d.txt");
intgral_1d_nonuni_cpp = importdata("..\resource\test_result\integral\integral_1d_non_uniform.txt");
intgral_2d_cpp = importdata("..\resource\test_result\integral\integral_matrix.txt");
intgral_2d_nonuni_cpp = importdata("..\resource\test_result\integral\integral_matrix_non_uniform.txt");