clc; clear;

out_dir = fullfile('..', 'test_result', 'matlab', 'interp');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

x = [0.0, 1.0, 2.0, 3.0];
y = [1.0, 3.0, 5.0, 7.0];
xq = [-1.0, 0.5, 1.5, 3.0, 4.0];
near_xq = [0.25, 1.4, 2.6];

linear = interp1(x, y, xq, 'linear', 'extrap');
spline = interp1(x, y, xq, 'spline', 'extrap');
pchip = interp1(x, y, xq, 'pchip', 'extrap');
makima = interp1(x, y, xq, 'makima', 'extrap');
nearest = interp1(x, y, near_xq, 'nearest', 'extrap');

writematrix([xq(:), linear(:)], ...
    fullfile(out_dir, 'linear.txt'), 'Delimiter', ' ');
writematrix([xq(:), spline(:)], ...
    fullfile(out_dir, 'spline.txt'), 'Delimiter', ' ');
writematrix([xq(:), pchip(:)], ...
    fullfile(out_dir, 'pchip.txt'), 'Delimiter', ' ');
writematrix([xq(:), makima(:)], ...
    fullfile(out_dir, 'makima.txt'), 'Delimiter', ' ');
writematrix([near_xq(:), nearest(:)], ...
    fullfile(out_dir, 'nearest.txt'), 'Delimiter', ' ');
