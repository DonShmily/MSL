script_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(script_dir, '..', '..', '..');
data_dir = fullfile(repo_root, 'test_result', 'difference', 'matlab_compare');

required = {'difference_vector.txt', 'difference_matrix.txt', ...
    'difference_row_diff.txt', 'difference_col_diff.txt'};
for i = 1:numel(required)
    if ~isfile(fullfile(data_dir, required{i}))
        error('Missing comparison data. Run "xmake run test_difference" first.');
    end
end

vec_data = readmatrix(fullfile(data_dir, 'difference_vector.txt'));
y = vec_data(:, 1);
diff_msl = vec_data(1:end-1, 2);
central_msl = vec_data(:, 3);

diff_ref = diff(y);
central_ref = [
    y(2) - y(1);
    (y(3:end) - y(1:end-2)) / 2.0;
    y(end) - y(end-1)
];

assert(max(abs(diff_msl - diff_ref)) < 1e-12, ...
    'Vector diff comparison failed.');
assert(max(abs(central_msl - central_ref)) < 1e-12, ...
    'Central gradient comparison failed.');

M = readmatrix(fullfile(data_dir, 'difference_matrix.txt'));
row_diff = readmatrix(fullfile(data_dir, 'difference_row_diff.txt'));
col_diff = readmatrix(fullfile(data_dir, 'difference_col_diff.txt'));

assert(max(abs(row_diff - diff(M, 1, 1)), [], 'all') < 1e-12, ...
    'Matrix row diff comparison failed.');
assert(max(abs(col_diff - diff(M, 1, 2)), [], 'all') < 1e-12, ...
    'Matrix col diff comparison failed.');

fprintf('Difference validation passed. vector diff error = %.3e\n', ...
    max(abs(diff_msl - diff_ref)));
