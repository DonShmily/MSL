clc;clear;close all
% compare_with_msl.m
% Run after C++ generated msl_out/*.csv
outdir = "msl_out";
if ~exist(outdir, 'dir')
    error('Directory "%s" not found. Run the C++ generator first.', outdir);
end

cases = {
    struct('tag','case01','type','low','N',2,'Wn',0.10),
    struct('tag','case02','type','low','N',4,'Wn',0.20),
    struct('tag','case03','type','low','N',6,'Wn',0.40),
    struct('tag','case04','type','low','N',4,'Wn',0.01),
    struct('tag','case05','type','high','N',2,'Wn',0.10),
    struct('tag','case06','type','high','N',4,'Wn',0.40),
    struct('tag','case07','type','band','N',2,'Wn',[0.10 0.15]),
    struct('tag','case08','type','band','N',3,'Wn',[0.05 0.40]),
    struct('tag','case09','type','stop','N',4,'Wn',[0.20 0.25]),
    struct('tag','case10','type','low','N',8,'Wn',0.48)
};

npoints = 1024;
freqs = linspace(0,1,npoints); % normalized 0..1 (1 => Nyquist, ω=π)

fprintf('Comparing MSL outputs in "%s" with MATLAB butter results\n\n', outdir);

for k = 1:length(cases)
    cs = cases{k};
    tag = cs.tag;
    fprintf('--- %s: %s N=%d Wn=%s ---\n', tag, cs.type, cs.N, mat2str(cs.Wn));

    % read MSL b,a
    b_msl = [];
    a_msl = [];
    try
        b_msl = csvread(fullfile(outdir, [tag '_b.csv']))';
        a_msl = csvread(fullfile(outdir, [tag '_a.csv']))';
    catch
        fprintf('  ERROR: failed to read %s files. Skipping.\n', tag);
        continue;
    end
    % ensure vectors
    b_msl = b_msl(:)'; a_msl = a_msl(:)';

    % compute matlab butter
    switch cs.type
        case 'low'
            [b_mat,a_mat] = butter(cs.N, cs.Wn, 'low');
        case 'high'
            [b_mat,a_mat] = butter(cs.N, cs.Wn, 'high');
        case 'band'
            [b_mat,a_mat] = butter(cs.N, cs.Wn, 'bandpass');
        case 'stop'
            [b_mat,a_mat] = butter(cs.N, cs.Wn, 'stop');
        otherwise
            error('unknown type %s', cs.type);
    end

    b_mat = b_mat(:)'; a_mat = a_mat(:)';

    % pad shorter vector with zeros for direct element comparison
    Lb = max(length(b_mat), length(b_msl));
    La = max(length(a_mat), length(a_msl));
    b_mat_p = [b_mat, zeros(1, Lb-length(b_mat))];
    b_msl_p = [b_msl, zeros(1, Lb-length(b_msl))];
    a_mat_p = [a_mat, zeros(1, La-length(a_mat))];
    a_msl_p = [a_msl, zeros(1, La-length(a_msl))];

    % numeric comparison
    diff_b = abs(b_mat_p - b_msl_p);
    diff_a = abs(a_mat_p - a_msl_p);
    maxdiff_b = max(diff_b);
    maxdiff_a = max(diff_a);
    rms_b = sqrt(mean(diff_b.^2));
    rms_a = sqrt(mean(diff_a.^2));

    % print small table
    fprintf('  matlab b (%d):\n   ', length(b_mat)); fprintf('%g ', b_mat); fprintf('\n');
    fprintf('  msl    b (%d):\n   ', length(b_msl)); fprintf('%g ', b_msl); fprintf('\n');
    fprintf('  diff b (padded to %d): max=%g, rms=%g\n', Lb, maxdiff_b, rms_b);

    fprintf('  matlab a (%d):\n   ', length(a_mat)); fprintf('%g ', a_mat); fprintf('\n');
    fprintf('  msl    a (%d):\n   ', length(a_msl)); fprintf('%g ', a_msl); fprintf('\n');
    fprintf('  diff a (padded to %d): max=%g, rms=%g\n', La, maxdiff_a, rms_a);

    % frequency response comparison
    H_mat = freqz(b_mat, a_mat, npoints);
    H_msl = freqz(b_msl, a_msl, npoints); % freqz accepts vectors of different lengths if a shorter -> zero pad

    mag_mat = abs(H_mat);
    mag_msl = abs(H_msl);
    mag_diff = abs(mag_mat - mag_msl);
    max_mag_diff = max(mag_diff);
    rms_mag = sqrt(mean(mag_diff.^2));

    fprintf('  freq response: max_mag_diff=%g, rms_mag=%g\n\n', max_mag_diff, rms_mag);

    % Plot coefficients and responses
    figure('Name',sprintf('%s coeffs and response',tag),'Visible','on');

    subplot(3,1,1);
    stem(0:Lb-1, b_mat_p, 'b','DisplayName','matlab b'); hold on;
    stem(0:Lb-1, b_msl_p, 'r','DisplayName','msl b','Marker','none');
    hold off; legend; title([tag ' numerator b (padded)']);

    subplot(3,1,2);
    stem(0:La-1, a_mat_p, 'b','DisplayName','matlab a'); hold on;
    stem(0:La-1, a_msl_p, 'r','DisplayName','msl a','Marker','none');
    hold off; legend; title([tag ' denominator a (padded)']);

    subplot(3,1,3);
    plot(freqs, 20*log10(mag_mat+eps), 'b', 'DisplayName','matlab |H| (dB)'); hold on;
    plot(freqs, 20*log10(mag_msl+eps), 'r--', 'DisplayName','msl |H| (dB)');
    plot(freqs, 20*log10(mag_diff+eps), 'k:', 'DisplayName','|Δ| (linear)'); % show raw diff in dB scale approx
    hold off; legend; xlabel('Normalized freq (0..1)'); ylabel('Magnitude (dB)'); title([tag ' freq response']);
end

fprintf('Done comparison. Inspect figures for details.\n');
