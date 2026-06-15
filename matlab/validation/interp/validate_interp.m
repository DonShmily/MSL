script_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(script_dir, '..', '..', '..');
data_dir = fullfile(repo_root, 'test_result', 'interp', 'matlab_compare');

interp_file = fullfile(data_dir, 'interp_results.txt');
nearest_file = fullfile(data_dir, 'interp_nearest.txt');
if ~isfile(interp_file) || ~isfile(nearest_file)
    error('Missing comparison data. Run "xmake run test_interp" first.');
end

x = [0.0; 1.0; 2.0; 3.0];
y = [1.0; 3.0; 5.0; 7.0];
data = readmatrix(interp_file);
xq = data(:, 1);

linear_ref = interp1(x, y, xq, 'linear', 'extrap');
spline_ref = interp1(x, y, xq, 'spline', 'extrap');
pchip_ref = interp1(x, y, xq, 'pchip', 'extrap');
makima_ref = interp1(x, y, xq, 'makima', 'extrap');

assert(max(abs(data(:, 2) - linear_ref)) < 1e-12, ...
    'Linear interpolation comparison failed.');
assert(max(abs(data(:, 3) - spline_ref)) < 1e-10, ...
    'Spline interpolation comparison failed.');
assert(max(abs(data(:, 4) - pchip_ref)) < 1e-12, ...
    'PCHIP interpolation comparison failed.');
assert(max(abs(data(:, 5) - makima_ref)) < 1e-12, ...
    'Akima/MAKIMA interpolation comparison failed.');

nearest = readmatrix(nearest_file);
nearest_ref = interp1(x, y, nearest(:, 1), 'nearest', 'extrap');
assert(max(abs(nearest(:, 2) - nearest_ref)) < 1e-12, ...
    'Nearest interpolation comparison failed.');

fprintf('Interp validation passed. linear error = %.3e\n', ...
    max(abs(data(:, 2) - linear_ref)));
