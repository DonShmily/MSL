clc; clear;

out_dir = fullfile('..', 'test_result', 'matlab', 'integral');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

y = [0.0; 1.0; 4.0; 9.0];
x = [0.0; 0.5; 2.0; 3.0];
parabola = [0.0; 1.0; 4.0];

cumtrapz_uniform = cumtrapz(y);
trapz_uniform = trapz(y);
cumtrapz_nonuniform = cumtrapz(x, y);
trapz_nonuniform = trapz(x, y);
simpson_three_points = (parabola(1) + 4.0 * parabola(2) + parabola(3)) / 3.0;

mat = [(0:3)', 2.0 * (0:3)'];
matrix_trapz = trapz(mat);
matrix_cumtrapz = cumtrapz(mat);

writematrix(cumtrapz_uniform, ...
    fullfile(out_dir, 'cumtrapz_uniform.txt'), 'Delimiter', ' ');
writematrix(trapz_uniform, ...
    fullfile(out_dir, 'trapz_uniform.txt'), 'Delimiter', ' ');
writematrix(cumtrapz_nonuniform, ...
    fullfile(out_dir, 'cumtrapz_nonuniform.txt'), 'Delimiter', ' ');
writematrix(trapz_nonuniform, ...
    fullfile(out_dir, 'trapz_nonuniform.txt'), 'Delimiter', ' ');
writematrix(simpson_three_points, ...
    fullfile(out_dir, 'simpson_three_points.txt'), 'Delimiter', ' ');
writematrix(matrix_trapz, ...
    fullfile(out_dir, 'matrix_trapz.txt'), 'Delimiter', ' ');
writematrix(matrix_cumtrapz, ...
    fullfile(out_dir, 'matrix_cumtrapz.txt'), 'Delimiter', ' ');
