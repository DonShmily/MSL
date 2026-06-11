clc; clear;

out_dir = fullfile('..', 'test_result', 'matlab', 'polynomial');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

coeffs_ascending = [1.0, 2.0, 3.0];
coeffs_descending = fliplr(coeffs_ascending);
query = [0.0, 1.0, 2.0];
values = polyval(coeffs_descending, query);
derivative_at_2 = polyval(polyder(coeffs_descending), 2.0);

x = [-2.0, -1.0, 0.0, 1.0, 2.0];
y = 1.0 - 2.0 .* x + 0.5 .* x .* x;
fit_descending = polyfit(x, y, 2);
fit_ascending = fliplr(fit_descending);
constant_fit = mean(y);

writematrix(values(:), fullfile(out_dir, 'polyval.txt'), 'Delimiter', ' ');
writematrix(derivative_at_2, ...
    fullfile(out_dir, 'derivative_at_2.txt'), 'Delimiter', ' ');
writematrix(fit_ascending(:), ...
    fullfile(out_dir, 'polyfit_coeffs_ascending.txt'), 'Delimiter', ' ');
writematrix(constant_fit, ...
    fullfile(out_dir, 'constant_fit.txt'), 'Delimiter', ' ');
