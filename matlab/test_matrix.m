clc;clear;
% 2025.11.08测试通过

%% 测试实矩阵SVD
A = [3,4;2,5;1,6];
[UA,SA,VA] = svd(A);
disp('A = '); disp(A);
disp('U = '); disp(UA);
disp('S = '); disp(SA);
disp('V = '); disp(VA);
disp('Reconstructed A = '); disp(UA*SA*VA');

%% 测试复矩阵SVD
B = [3+1i,4-1i;2-1i,5+1i;1+1i,6-1i];
[UB,SB,VB] = svd(B);
disp('B = '); disp(B);
disp('U = '); disp(UB);
disp('S = '); disp(SB);
disp('V = '); disp(VB);
disp('Reconstructed B = '); disp(UB*SB*VB');