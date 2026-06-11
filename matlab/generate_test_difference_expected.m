clc; clear;

out_dir = fullfile('..', 'test_result', 'matlab', 'difference');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

y = [1.0; 4.0; 9.0; 16.0];
diff_1d = diff(y);
forward_gradient_dx_0_5 = diff(y) ./ 0.5;
central_gradient = gradient(y, 1.0);

mat = zeros(3, 3);
for i = 1:3
    for j = 1:3
        mat(i, j) = 10.0 * (i - 1) + (j - 1);
    end
end
row_diff = diff(mat, 1, 1);
col_diff = diff(mat, 1, 2);
laplacian_center = 0.0;

writematrix(diff_1d, fullfile(out_dir, 'diff_1d.txt'), 'Delimiter', ' ');
writematrix(forward_gradient_dx_0_5, ...
    fullfile(out_dir, 'forward_gradient_dx_0_5.txt'), 'Delimiter', ' ');
writematrix(central_gradient, ...
    fullfile(out_dir, 'central_gradient.txt'), 'Delimiter', ' ');
writematrix(row_diff, fullfile(out_dir, 'row_diff.txt'), 'Delimiter', ' ');
writematrix(col_diff, fullfile(out_dir, 'col_diff.txt'), 'Delimiter', ' ');
writematrix(laplacian_center, ...
    fullfile(out_dir, 'laplacian_center.txt'), 'Delimiter', ' ');
