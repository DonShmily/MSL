clc; clear;

out_dir = fullfile('..', 'test_result', 'matlab', 'matrix');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

A = [3.0, 4.0; 2.0, 5.0; 1.0, 6.0];
[U, S, V] = svd(A, 'econ');
reconstructed = U * S * V';

B = [1.0, 11.0, 21.0; 2.0, 12.0, 22.0];
C = [1.0, 11.0; 2.0, 12.0; 3.0, 13.0];
matmul = B * C;

writematrix(A, fullfile(out_dir, 'svd_input.txt'), 'Delimiter', ' ');
writematrix(U, fullfile(out_dir, 'svd_u.txt'), 'Delimiter', ' ');
writematrix(S, fullfile(out_dir, 'svd_s.txt'), 'Delimiter', ' ');
writematrix(V', fullfile(out_dir, 'svd_vt.txt'), 'Delimiter', ' ');
writematrix(reconstructed, ...
    fullfile(out_dir, 'svd_reconstructed.txt'), 'Delimiter', ' ');
writematrix(matmul, fullfile(out_dir, 'matmul.txt'), 'Delimiter', ' ');
