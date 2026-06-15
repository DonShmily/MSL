script_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(script_dir, '..', '..', '..');
data_dir = fullfile(repo_root, 'test_result', 'integral', 'matlab_compare');

summary_file = fullfile(data_dir, 'integral_summary.txt');
cum_file = fullfile(data_dir, 'integral_cumulative.txt');
if ~isfile(summary_file) || ~isfile(cum_file)
    error('Missing comparison data. Run "xmake run test_integral" first.');
end

summary = readmatrix(summary_file);
cum_data = readmatrix(cum_file);

x = cum_data(:, 1);
y = cum_data(:, 2);
cum_uniform = cum_data(:, 3);
cum_nonuniform = cum_data(:, 4);

expected = [
    trapz(y);
    trapz(x, y);
    integral(@(xx) xx.^2, 0.0, 2.0);
    integral(@(xx) xx.^2, 0.0, 2.0);
    integral(@(xx) xx.^2, 0.0, 4.0);
    integral(@sin, 0.0, pi)
];

assert(max(abs(summary - expected)) < 1e-10, ...
    'Integral summary comparison failed.');
assert(max(abs(cum_uniform - cumtrapz(y))) < 1e-12, ...
    'Uniform cumtrapz comparison failed.');
assert(max(abs(cum_nonuniform - cumtrapz(x, y))) < 1e-12, ...
    'Non-uniform cumtrapz comparison failed.');

fprintf('Integral validation passed. summary error = %.3e\n', ...
    max(abs(summary - expected)));
