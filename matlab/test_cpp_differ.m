clc;clear;close all
% 2025.10.23测试通过

ori_data = importdata("..\resource\KunmingSSJY.txt");

diff_1d = [0;diff(ori_data(:,1))];
diff_2d_row = [zeros(3e4,1), diff(ori_data,1,2)];
diff_2d_col = [zeros(1,6);diff(ori_data,1,1)];
grad_1d_row = gradient(ori_data(:,1));
[grad1, grad2] = gradient(ori_data);

diff_1d_cpp = importdata("..\resource\test_result\difference\diff_1d.txt");
diff_2d_row_cpp = importdata("..\resource\test_result\difference\diff_2d_row.txt");
diff_2d_col_cpp = importdata("..\resource\test_result\difference\diff_2d_col.txt");
grad_1d_row_cpp = importdata("..\resource\test_result\difference\cent_grad_1d.txt");
grad_2d_row_cpp = importdata("..\resource\test_result\difference\cent_grad_2d_row.txt");
grad_2d_x_cpp = importdata("..\resource\test_result\difference\cent_grad_2d_col_x.txt");
grad_2d_y_cpp = importdata("..\resource\test_result\difference\cent_grad_2d_col_y.txt");