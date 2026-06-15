script_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(script_dir, '..', '..', '..');
data_dir = fullfile(repo_root, 'test_result', 'matrix', 'matlab_compare');

required = {'matrix_A.txt', 'matrix_B.txt', 'matrix_product.txt', ...
    'matrix_svd_input.txt', 'matrix_svd_reconstruction.txt'};
for i = 1:numel(required)
    if ~isfile(fullfile(data_dir, required{i}))
        error('Missing comparison data. Run "xmake run test_matrix" first.');
    end
end

A = readmatrix(fullfile(data_dir, 'matrix_A.txt'));
B = readmatrix(fullfile(data_dir, 'matrix_B.txt'));
C = readmatrix(fullfile(data_dir, 'matrix_product.txt'));
SVD_input = readmatrix(fullfile(data_dir, 'matrix_svd_input.txt'));
SVD_reconstructed = readmatrix(fullfile(data_dir, 'matrix_svd_reconstruction.txt'));

assert(max(abs(C - A * B), [], 'all') < 1e-12, ...
    'Matrix product comparison failed.');
assert(max(abs(SVD_reconstructed - SVD_input), [], 'all') < 1e-10, ...
    'Matrix SVD reconstruction comparison failed.');

fprintf('Matrix validation passed. product error = %.3e\n', ...
    max(abs(C - A * B), [], 'all'));
