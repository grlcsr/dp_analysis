%% Parameters
clear; clc; close all;

figure_idx = 0;
plot_figs = 1;

entropies  = {"M4", "M8", "M16", "M32", "M64"};

bits           = {"clear", "set", "flip", "alternate", "correlation"};
bits_ltx       = {"CLEAR", "SET", "FLIP", "ALTERNATE", "CORRELATION"};

bin_names          = {'rap2', 'openssl', 'xor', 'mgb_bad'};
libraries          = {'opendp', 'opendp_mod', 'ibm'};

iterations     = 1000000;

bin            = bin_names{1};
library        = libraries{3};

epsilon        = 0.1;
bin_sizes      = [1,5];
hs             = 1:length(entropies);
bs             = 1:length(bits);

dbsize1        = 10000;
dbsize2        = 10001;
midpoint = (dbsize1 + dbsize2)/2;

%% First study the plots with no entropy manipulation
mat_name = "./data_mats/eps_" + sprintf('%.2f', epsilon) + "/" + library + "/" + bin + "/0H.mat";
mat      = load(mat_name);
[x1, x2] = load_counts(mat, 1, iterations);

for bin_size = bin_sizes
    [min_v, max_v, edges, bin_c] = make_edges_from_data(x1, x2, bin_size, midpoint);

    % Histograms + binomial errors
    [bin_v1, p_err1] = hist_with_errs(x1, edges);
    [bin_v2, p_err2] = hist_with_errs(x2, edges);

    % Distributions
    if plot_figs > 0
        [figure_idx] = plot_distributions(bin_c, bin_v1/iterations, p_err1/iterations, bin_v2/iterations, p_err2/iterations, ...
            ["DP Definition - Distributions", "No H manipulation - $\epsilon$ = " + epsilon + " - N = " + iterations + " - bin size = " + bin_size], ...
            figure_idx);
    end

    % Ratios + theory
    [ratio_bin, ratio_bin_err, bin_c] = ratio_and_err(bin_v1, bin_v2, bin_c);
    r_theory = exp(epsilon * sign(midpoint - bin_c));

    if plot_figs > 0
        [figure_idx] = plot_privacy_loss(bin_c, ratio_bin, ratio_bin_err, r_theory, ...
            ["DP Definition - Privacy Loss rv.", "No H manipulation - $\epsilon$ = " + epsilon + " - N = " + iterations + " - bin size = " + bin_size], ...
            figure_idx, epsilon, [bin_c(1), bin_c(end)]);
    end

    fprintf("\nAnalysis on H = 0 (no entropy manipulation) on %s samples with bin size %d\n", library, bin_size);
    chi2_report(bin_c, ratio_bin, ratio_bin_err, r_theory);
    chi_report_68(bin_c, ratio_bin, ratio_bin_err, r_theory);

    % Sign test
    sign_test(bin_c, ratio_bin, ratio_bin_err, midpoint, epsilon);
    sign_test_68(bin_c, ratio_bin, ratio_bin_err, midpoint, epsilon);
end

%% Manipulated entropy loops
for h = hs
    cur_h = entropies{h};
    for b = bs
        cur_b = bits{b};

        mat_name = "./data_mats/eps_" + sprintf('%.2f', epsilon) + "/" + library + "/" + bin + "/" + cur_h + "_" + cur_b + ".mat";
        mat      = load(mat_name);
        [x1, x2] = load_counts(mat, 1, iterations);

        for bin_size = bin_sizes
            [min_v, max_v, edges, bin_c] = make_edges_from_data(x1, x2, bin_size, midpoint);

            % Histograms + binomial errors
            [bin_v1, p_err1] = hist_with_errs(x1, edges);
            [bin_v2, p_err2] = hist_with_errs(x2, edges);

            % Distributions
            if plot_figs > 0
                [figure_idx] = plot_distributions(bin_c, bin_v1/iterations, p_err1/iterations, bin_v2/iterations, p_err2/iterations, ...
                    ["DP Definition - Privacy Loss rv.", entropies{h} + " manipulation " + bits_ltx{b} + " - $\epsilon$ = " + epsilon + " - N = " + iterations + " - bin size = " + bin_size], ...
                    figure_idx);
            end

            % Ratios + theory
            [ratio_bin, ratio_bin_err, bin_c] = ratio_and_err(bin_v1, bin_v2, bin_c);
            midpoint = (dbsize1 + dbsize2)/2;
            r_theory = exp(epsilon * sign(midpoint - bin_c));

            if plot_figs > 0
                [figure_idx] = plot_privacy_loss(bin_c, ratio_bin, ratio_bin_err, r_theory, ...
                    ["DP Definition - Privacy Loss rv.", entropies{h} + " manipulation " + bits_ltx{b} + " - $\epsilon$ = " + epsilon + " - N = " + iterations + " - bin size = " + bin_size], ...
                    figure_idx, epsilon, [bin_c(1), bin_c(end)]);
            end

            fprintf("\nAnalysis on \'%s %s\' entropy manipulation on %s samples with bin size %d\n", cur_h, cur_b, library, bin_size);
            % Fit with function fixed at theoretical epsilon and compute χ²
            chi2_report(bin_c, ratio_bin, ratio_bin_err, r_theory);
            chi_report_68(bin_c, ratio_bin, ratio_bin_err, r_theory);

            % Sign test
            sign_test(bin_c, ratio_bin, ratio_bin_err, midpoint, epsilon)
            sign_test_68(bin_c, ratio_bin, ratio_bin_err, midpoint, epsilon)
        end
    end
end

%% Functions
function [x1, x2] = load_counts(mat, start, fin)
x1 = mat.dp1(start:fin);
x2 = mat.dp2(start:fin);
end

function [min_v, max_v, edges, centers] = make_edges_from_data(x1, x2, step, midpoint)
min_v   = min(min(x1, x2));
max_v   = max(max(x1, x2));

% we want an edge to be exactly on the step point to avoid skewing the
% samples
c1      = midpoint - step/2;
cs      = sort(c1:-step:min_v);
centers = [cs, c1 + step:step:max_v];
edges   = centers(1) - step/2 : step : centers(end) + step/2;
end

function [bin_v, p_err] = hist_with_errs(x, edges)
bin_v = histcounts(x, edges);
[~, p_err] = Binomial_err(bin_v);
end

function [ratio_v, ratio_err, bin_c] = ratio_and_err(v1, v2, bin_c)
assert(numel(v1)==numel(v2) && numel(v1)==numel(bin_c), ...
    'v1, v2, and bin_c must have the same length.');

% Only compute where both bins are nonzero and finite
nz = (v1 ~= 0) & (v2 ~= 0) & isfinite(v1) & isfinite(v2);

r = nan(size(v1));
e = nan(size(v1));

v1n = v1(nz);
v2n = v2(nz);

r(nz) = v1n ./ v2n;
e(nz) = sqrt( v1n ./ (v2n.^2) + (v1n.^2) ./ (v2n.^3) );

% Keep only finite results and align bin_c
valid = isfinite(r) & isfinite(e);

ratio_v   = r(valid);
ratio_err = e(valid);
bin_c     = bin_c(valid);

ratio_v   = ratio_v(:).';
ratio_err = ratio_err(:).';
bin_c     = bin_c(:).';
end

function figure_idx = plot_distributions(bin_c, v1, e1, v2, e2, title_lines, figure_idx)
figure_idx = figure_idx + 1;
figure(figure_idx); hold on; grid on; box on;
errorbar(bin_c, v1, e1, '*');
errorbar(bin_c, v2, e2, '*');
title(title_lines, 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Noise'); ylabel('Samples');
end

function figure_idx = plot_privacy_loss(bin_c, ratio_v, ratio_err, r_theory, title_lines, figure_idx, epsilon, xlim_override)
if nargin < 7, xlim_override = [bin_c(1), bin_c(end)]; end
figure_idx = figure_idx + 1;
figure(figure_idx); hold on; grid on; box on;
errorbar(bin_c, ratio_v, ratio_err, '*');
plot(bin_c, r_theory, '--', 'LineWidth', 1.2);
yline(1,  'LineStyle', '--');
yline(exp(epsilon),   '-.', 'Color', [1, 0.5, 0.3125]);
yline(exp(-epsilon),  '-.', 'Color', [1, 0.5, 0.3125]);
yline(1, 'LineStyle', '--');
title(title_lines, 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Privacy Loss'); xlabel('Samples');
xlim(xlim_override);
end

function chi2_report(bin_c, ratio_bin, ratio_bin_err, r_theory)
ndf = 0;
chi2_fit = 0;

% ---------- chi2 in spazio lineare ----------
for j = 1:length(bin_c)
    if (ratio_bin_err(j) ~= 0 && ...
        ~isnan(ratio_bin_err(j)) && ...
        ~isnan(ratio_bin(j)) && ...
        ~isnan(r_theory(j)))

        chi2_fit = chi2_fit + ((r_theory(j) - ratio_bin(j)) / ratio_bin_err(j))^2;
        ndf = ndf + 1;
    end
end

chi2_red     = chi2_fit / ndf;
chi2_red_std = sqrt(2 / ndf);
pval         = 1 - chi2cdf(chi2_fit, ndf);

fprintf("  χ² fit with theoretical epsilon (linear space):\n");
fprintf("  χ²_red = %.5f ± %.5f   (ndf = %d), p-value = %.3g\n", ...
        chi2_red, chi2_red_std, ndf, pval);

% ---------- chi2 in log-spazio ----------
chi2_fit_log = 0;
ndf_log      = 0;

for j = 1:length(bin_c)
    % condizioni: errori definiti, ratio > 0, teoria > 0
    if (ratio_bin_err(j) ~= 0 && ...
        ~isnan(ratio_bin_err(j)) && ...
        ~isnan(ratio_bin(j)) && ...
        ~isnan(r_theory(j)) && ...
        ratio_bin(j) > 0 && ...
        r_theory(j) > 0)

        % errore in log-spazio: σ_log ≈ σ_r / r
        sigma_log = ratio_bin_err(j) / ratio_bin(j);

        if sigma_log > 0 && ~isnan(sigma_log) && ~isinf(sigma_log)
            diff_log    = log(r_theory(j)) - log(ratio_bin(j));
            chi2_fit_log = chi2_fit_log + (diff_log / sigma_log)^2;
            ndf_log      = ndf_log + 1;
        end
    end
end

chi2_red_log     = chi2_fit_log / ndf_log;
chi2_red_log_std = sqrt(2 / ndf_log);
pval_log         = 1 - chi2cdf(chi2_fit_log, ndf_log);

fprintf("  χ² fit with theoretical epsilon (log-space):\n");
fprintf("  χ²_red = %.5f ± %.5f   (ndf = %d), p-value = %.3g\n", ...
        chi2_red_log, chi2_red_log_std, ndf_log, pval_log);
end


function sign_test(bin_c, ratio_bin, ratio_bin_err, midpoint, epsilon)
% Sign test - both wrt r_theory and sample per sample
% Matlab function
idx_lower = bin_c < midpoint;
s = ones(size(bin_c)); s(~idx_lower) = -1;

dlog = log(ratio_bin) - s .* epsilon;
valid = isfinite(dlog);
dlog = dlog(valid);
pval = signtest(dlog, 0);
% fprintf('  Sign test (two tails, log-space): n=%d, p-value=%.4g\n', numel(dlog), pval);

% Option one: wrt r_theory not considering error bars
d = zeros(size(ratio_bin));
d(idx_lower) = ratio_bin(idx_lower) - exp(epsilon);
d(~idx_lower) = exp(-epsilon) - ratio_bin(~idx_lower);

is_tie = ~isfinite(d) | (d == 0); % need to remove zeros (and infinites/NaNs in case)
signs  = sign(d(~is_tie));

k_pos  = sum(signs > 0);
n  = numel(signs);
k_min  = min(k_pos, n - k_pos);

pval  = min(2 * binocdf(k_min, n, 0.5), 1);
% fprintf('  Sign test (two tails): n = %d, k_pos = %d, p-value = %.4g \n', n, k_pos, pval);

% Option 2: we consider error bars -> not really a sign test cuz we
% compare how many values are within the ksigma range with expected
% binomial
for k = 1:2
    t = nan(size(ratio_bin));
    t(idx_lower)  = exp(epsilon);
    t(~idx_lower) = exp(-epsilon);

    z = (ratio_bin - t) ./ ratio_bin_err;

    valid = isfinite(z) & isfinite(ratio_bin_err) & (ratio_bin_err > 0);
    z = z(valid);
    % z(~idx_lower) = z(~idx_lower)*-1;

    out = abs(z) > k;
    % out = z > k;
    k_out = sum(out);
    n = numel(z);
    p0 = 2 * (1 - normcdf(k));
    pval  = binocdf(k_out, n, p0, "upper");

    fprintf('  %gσ exceedance test: k_out = %d/%d, p-value = %.3g\n', k, k_out, n, pval);

    phat = k_out / n;
    [~, ci] = binofit(k_out, n);   % 95% CI for exceedance rate
    fprintf('  Observed rate = %.3f (95%% CI [%.3f, %.3f]); expected under H0 = %.4f\n', ...
        phat, ci(1), ci(2), p0);

    % Running the same test in log space

    s = ones(size(ratio_bin)); s(~idx_lower) = -1;
    dlog = log(ratio_bin) - s.*epsilon;

    % delta-method SE for log(r): sigma_log ≈ sigma_ratio / ratio
    sigma_log = ratio_bin_err ./ ratio_bin;

    z_log = dlog ./ sigma_log;
    valid = isfinite(z_log) & isfinite(sigma_log) & (sigma_log > 0);
    z_log = z_log(valid);

    out_log   = abs(z_log) > k;
    k_out_log = sum(out_log);
    n_log     = numel(z_log);
    p0        = 2 * (1 - normcdf(k));
    pval_log  = binocdf(k_out_log, n_log, p0, 'upper');

    fprintf('  %gσ exceedance (log-space): k_out = %d/%d, p-value = %.3g\n', ...
        k, k_out_log, n_log, pval_log);
end
end

function chi_report_68(bin_c, ratio_bin, ratio_bin_err, r_theory)

N = numel(bin_c);

% Number of elements to keep (central 68%)
K = max(1, round(0.68 * N));

% index of element closest to median
[~, center_idx] = min(abs(bin_c - median(bin_c)));

% compute start and stop indices trying to balance left/right
half_left  = floor((K-1)/2);
half_right = K-1-half_left;

start_idx = center_idx - half_left;
end_idx   = center_idx + half_right;

if start_idx < 1
    shift = 1 - start_idx;
    start_idx = 1;
    end_idx = min(N, end_idx + shift);
end
if end_idx > N
    shift = end_idx - N;
    end_idx = N;
    start_idx = max(1, start_idx - shift);
end

sel_bin_c        = bin_c(start_idx:end_idx);
sel_ratio_bin    = ratio_bin(start_idx:end_idx);
sel_ratio_bin_err = ratio_bin_err(start_idx:end_idx);
sel_r_theory = r_theory(start_idx:end_idx);

fprintf("χ²_red fit for the central 68%%\n");
chi2_report(sel_bin_c, sel_ratio_bin, sel_ratio_bin_err, sel_r_theory);
end

function sign_test_68(bin_c, ratio_bin, ratio_bin_err, midpoint, epsilon)
    
    N = numel(bin_c);
    
    % Number of elements to keep (central 68%)
    K = max(1, round(0.68 * N));
    
    % index of element closest to median
    [~, center_idx] = min(abs(bin_c - midpoint));
    
    % compute start and stop indices trying to balance left/right
    half_left  = floor((K-1)/2);
    half_right = K-1-half_left;
    
    start_idx = center_idx - half_left;
    end_idx   = center_idx + half_right;
    
    if start_idx < 1
        shift = 1 - start_idx;
        start_idx = 1;
        end_idx = min(N, end_idx + shift);
    end
    if end_idx > N
        shift = end_idx - N;
        end_idx = N;
        start_idx = max(1, start_idx - shift);
    end
    
    sel_bin_c        = bin_c(start_idx:end_idx);
    sel_ratio_bin    = ratio_bin(start_idx:end_idx);
    sel_ratio_bin_err = ratio_bin_err(start_idx:end_idx);
    
    fprintf("sigma exceedance test for the central 68%%\n");
    sign_test(sel_bin_c, sel_ratio_bin, sel_ratio_bin_err, midpoint, epsilon)

end