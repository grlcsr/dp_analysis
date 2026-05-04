%% Parameters
clear; clc; close all;

figure_idx = 0;
plot_figs = 1;
jump_bins = 1; % if set, exclude from tests the two bins at center of distributions (jump of the plrv)
print_68 = 0; % if set, also print analysis for the central 68% of bins

entropies = {"M4", "M8", "M16", "M32", "M64"};

bits           = {"bit_flip", "alternate_set", "correlation"};
bits_ltx       = {"FLIP", "ALTERNATE", "CORRELATION"};

bin_names          = {'rap2', 'openssl', 'mgb_bad'};
libraries          = {'ibm'};

iterations     = 1000000;

bin            = bin_names{1};
library        = libraries{3};

epsilon        = 0.1;
bin_sizes      = [1];%, 3, 5];
hs             = 1:length(entropies);
bs             = 1:length(bits);

dbsize1        = 10000;
dbsize2        = 10001;
midpoint = (dbsize1 + dbsize2)/2;

n_bins_side = 75; % 75 bins left + 75 bins right of midpoint

bin

%% First study the plots with no entropy manipulation
mat_name = "./data_mats/eps_" + sprintf('%.2f', epsilon) + "/" + library + "/" + bin + "/0H.mat";
mat      = load(mat_name);
[x1, x2] = load_counts(mat, 1, iterations);

for bin_size = bin_sizes
    [edges, bin_c] = make_fixed_edges(midpoint, bin_size, n_bins_side);

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
    [ratio_bin, ratio_bin_err, bin_c_valid] = ratio_and_err(bin_v1, bin_v2, bin_c);
    r_theory = exp(epsilon * sign(midpoint - bin_c_valid));

    if plot_figs > 0
        [figure_idx] = plot_privacy_loss(bin_c_valid, ratio_bin, ratio_bin_err, r_theory, ...
            ["DP Definition - Privacy Loss rv.", "No H manipulation - $\epsilon$ = " + epsilon + " - N = " + iterations + " - bin size = " + bin_size], ...
            figure_idx, epsilon, [bin_c_valid(1), bin_c_valid(end)]);
    end

    fprintf("\nAnalysis on H = 0 (no entropy manipulation) on %s samples with bin size %d\n", library, bin_size);
    run_all_tests(bin_c_valid, ratio_bin, ratio_bin_err, r_theory, midpoint, epsilon, dbsize1, dbsize2, jump_bins, print_68);
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
            [edges, bin_c] = make_fixed_edges(midpoint, bin_size, n_bins_side);

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
            [ratio_bin, ratio_bin_err, bin_c_valid] = ratio_and_err(bin_v1, bin_v2, bin_c);
            r_theory = exp(epsilon * sign(midpoint - bin_c_valid));

            if plot_figs > 0
                [figure_idx] = plot_privacy_loss(bin_c_valid, ratio_bin, ratio_bin_err, r_theory, ...
                    ["DP Definition - Privacy Loss rv.", entropies{h} + " manipulation " + bits_ltx{b} + " - $\epsilon$ = " + epsilon + " - N = " + iterations + " - bin size = " + bin_size], ...
                    figure_idx, epsilon, [bin_c_valid(1), bin_c_valid(end)]);
            end

            fprintf("\nAnalysis on \'%s %s\' entropy manipulation on %s samples with bin size %d\n", cur_h, cur_b, library, bin_size);
            run_all_tests(bin_c_valid, ratio_bin, ratio_bin_err, r_theory, midpoint, epsilon, dbsize1, dbsize2, jump_bins, print_68);
        end
    end
end

%% Functions

function run_all_tests(bin_c, ratio_bin, ratio_bin_err, r_theory, midpoint, epsilon, d1, d2, jump_bins, print_68)
% Run all tests on full bins
% Optionally, repeat after removing jump bins or run on central 68%.

    % Full histogram
    fprintf("  [Full histogram]\n");
    chi2_report(bin_c, ratio_bin, ratio_bin_err, r_theory);
    sign_test(bin_c, ratio_bin, ratio_bin_err, midpoint, epsilon);

    % Full histogram, jump bins excluded
    if jump_bins
        [bin_c_nj, ratio_nj, err_nj, th_nj] = remove_jump_bins(bin_c, ratio_bin, ratio_bin_err, r_theory, d1, d2);
        fprintf("  [Full histogram, jump bins excluded]\n");
        chi2_report(bin_c_nj, ratio_nj, err_nj, th_nj);
        sign_test(bin_c_nj, ratio_nj, err_nj, midpoint, epsilon);
    end

    if print_68
        % Central 68%
        [sel_c, sel_r, sel_e, sel_t] = select_central_68(bin_c, ratio_bin, ratio_bin_err, r_theory);
        fprintf("  [Central 68%%]\n");
        chi2_report(sel_c, sel_r, sel_e, sel_t);
        sign_test(sel_c, sel_r, sel_e, midpoint, epsilon);

        % Central 68%, jump bins excluded
        if jump_bins
            [sel_c_nj, sel_r_nj, sel_e_nj, sel_t_nj] = remove_jump_bins(sel_c, sel_r, sel_e, sel_t, d1, d2);
            fprintf("  [Central 68%%, jump bins excluded]\n");
            chi2_report(sel_c_nj, sel_r_nj, sel_e_nj, sel_t_nj);
            sign_test(sel_c_nj, sel_r_nj, sel_e_nj, midpoint, epsilon);
        end
    end
end

function [sel_c, sel_r, sel_e, sel_t] = select_central_68(bin_c, ratio_bin, ratio_bin_err, r_theory)
    N = numel(bin_c);
    K = max(1, round(0.68 * N));
    [~, center_idx] = min(abs(bin_c - median(bin_c)));

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

    sel_c = bin_c(start_idx:end_idx);
    sel_r = ratio_bin(start_idx:end_idx);
    sel_e = ratio_bin_err(start_idx:end_idx);
    sel_t = r_theory(start_idx:end_idx);
end

function [x1, x2] = load_counts(mat, start, fin)
    x1 = mat.dp1(start:fin);
    x2 = mat.dp2(start:fin);
end

function [edges, centers] = make_fixed_edges(midpoint, step, n_bins_side)
    left_centers  = midpoint - step/2 : -step : midpoint - (2*n_bins_side-1)*step/2;
    left_centers  = sort(left_centers);
    right_centers = midpoint + step/2 : step : midpoint + (2*n_bins_side-1)*step/2;
    centers = [left_centers, right_centers];
    edges   = (centers(1) - step/2) : step : (centers(end) + step/2);
end

function [bin_v, p_err] = hist_with_errs(x, edges)
    bin_v = histcounts(x, edges);
    [~, p_err] = Binomial_err(bin_v);
end

function [bin_c, ratio_v, ratio_err, r_theory] = remove_jump_bins(bin_c, ratio_v, ratio_err, r_theory, d1, d2)
    keep = (bin_c ~= d1) & (bin_c ~= d2);
    n_removed = sum(~keep);
    if n_removed > 0
        fprintf('    Removed %d jump bin(s) at centers:', n_removed);
        fprintf(' %.1f', bin_c(~keep));
        fprintf('\n');
    end
    bin_c     = bin_c(keep);
    ratio_v   = ratio_v(keep);
    ratio_err = ratio_err(keep);
    r_theory  = r_theory(keep);
end

function [ratio_v, ratio_err, bin_c] = ratio_and_err(v1, v2, bin_c)
    assert(numel(v1)==numel(v2) && numel(v1)==numel(bin_c), ...
        'v1, v2, and bin_c must have the same length.');

    n_input = numel(bin_c);

    nz = (v1 ~= 0) & (v2 ~= 0) & isfinite(v1) & isfinite(v2);

    r = nan(size(v1));
    e = nan(size(v1));

    v1n = v1(nz);
    v2n = v2(nz);

    r(nz) = v1n ./ v2n;
    e(nz) = sqrt( v1n ./ (v2n.^2) + (v1n.^2) ./ (v2n.^3) );

    valid = isfinite(r) & isfinite(e);

    ratio_v   = r(valid);
    ratio_err = e(valid);
    bin_c     = bin_c(valid);

    ratio_v   = ratio_v(:).';
    ratio_err = ratio_err(:).';
    bin_c     = bin_c(:).';

    n_dropped = n_input - numel(bin_c);
    if n_dropped > 0
        fprintf(2, '  WARNING: %d / %d bins dropped (zero counts or non-finite ratio). Surviving bins: %d\n', ...
            n_dropped, n_input, numel(bin_c));
    end
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
    title(title_lines, 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Privacy Loss'); xlabel('Samples');
    xlim(xlim_override);
end

function chi2_report(bin_c, ratio_bin, ratio_bin_err, r_theory)
    ndf = 0;
    chi2_fit = 0;

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

    fprintf("    χ²_red = %.5f ± %.5f   (ndf = %d), p-value = %.3g\n", ...
            chi2_red, chi2_red_std, ndf, pval);
end

function sign_test(bin_c, ratio_bin, ratio_bin_err, midpoint, epsilon)
    idx_lower = bin_c < midpoint;

    % Sign test
    d = zeros(size(ratio_bin));
    d(idx_lower) = ratio_bin(idx_lower) - exp(epsilon);
    d(~idx_lower) = exp(-epsilon) - ratio_bin(~idx_lower);

    is_tie = ~isfinite(d) | (d == 0);
    signs  = sign(d(~is_tie));

    k_pos  = sum(signs > 0);
    n  = numel(signs);
    k_min  = min(k_pos, n - k_pos);
    pval  = min(2 * binocdf(k_min, n, 0.5), 1);
    fprintf('    Sign test: n = %d, k_pos = %d, p-value = %.4g\n', n, k_pos, pval);

    % Sigma exceedance test
    for k = 1:2
        t = nan(size(ratio_bin));
        t(idx_lower)  = exp(epsilon);
        t(~idx_lower) = exp(-epsilon);

        z = (ratio_bin - t) ./ ratio_bin_err;

        valid = isfinite(z) & isfinite(ratio_bin_err) & (ratio_bin_err > 0);
        z = z(valid);

        out = abs(z) > k;
        k_out = sum(out);
        n = numel(z);
        p0 = 2 * (1 - normcdf(k));
        pval  = binocdf(k_out, n, p0, "upper");

        fprintf('    %gσ exceedance: k_out = %d/%d, p-value = %.3g\n', k, k_out, n, pval);

        phat = k_out / n;
        [~, ci] = binofit(k_out, n);
        fprintf('    Observed rate = %.3f (95%% CI [%.3f, %.3f]); expected = %.4f\n', ...
            phat, ci(1), ci(2), p0);
    end
end
