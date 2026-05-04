%% Parameters
clear; clc; close all;

figure_idx = 0;

libraries = {'ibm'};
bin_names = {'rap2', 'openssl', 'mgb_bad'};

iterations = 1000000;

bin     = bin_names{1};
library = libraries{2};

switch lower(library)
    case {'opendp', 'opendp_mod'}, tot_iterations = 100e6;
    case 'ibm',                    tot_iterations = 91e6;
    otherwise,                     tot_iterations = 0;
end
tot_data = floor(tot_iterations / iterations);

epsilon  = 0.1;
bin_size = 1;

dbsize1  = 10000;
dbsize2  = 10001;
midpoint = (dbsize1 + dbsize2) / 2;

mat_name = "./data_mats/eps_" + sprintf('%.2f', epsilon) + "/" ...
         + library + "/" + bin + "/0H.mat";
mat = load(mat_name);

%% First sample: histograms, ratio, plots
[x1, x2] = load_counts(mat, 1, iterations);

[~, ~, edges, bin_c] = make_edges_from_data(x1, x2, bin_size, midpoint);

[bin_v1, p_err1] = hist_with_errs(x1, edges);
[bin_v2, p_err2] = hist_with_errs(x2, edges);

[ratio_bin, ratio_bin_err] = ratio_and_err(bin_v1, bin_v2);
r_theory = exp(epsilon * sign(midpoint - bin_c));

figure_idx = plot_distributions(bin_c, bin_v1, p_err1, bin_v2, p_err2, ...
    ["DP Definition - Distributions", ...
     "No H manipulation - $\epsilon$ = " + epsilon ...
     + " - N = " + iterations + " - bin size = " + bin_size], ...
    figure_idx);

figure_idx = plot_privacy_loss(bin_c, ratio_bin, ratio_bin_err, r_theory, ...
    ["DP Definition - Privacy Loss rv.", ...
     "No H manipulation - $\epsilon$ = " + epsilon ...
     + " - N = " + iterations + " - bin size = " + bin_size], ...
    figure_idx, epsilon, [bin_c(1), bin_c(end)]);

%% Analytical error propagation
sigma_a_2 = var(x1) * ((-4*bin_v1.*bin_v2 + 2*bin_v1*iterations + bin_v2*iterations) ...
            ./ (bin_v2.^3 * iterations)).^2 ...
          + var(x2) * (bin_v1 .* (4*bin_v1.*bin_v2 - 3*bin_v1*iterations - 2*bin_v2*iterations) ...
            ./ (bin_v2.^4 * iterations)).^2;

%% Repeated samples for empirical standard deviation
n_samples = tot_data - 1;
sample_points_ratio     = zeros(n_samples, length(ratio_bin));
sample_points_ratio_err = zeros(n_samples, length(ratio_bin));

for j = 2:tot_data
    start = 1 + (j - 1) * iterations;
    fin   = start + iterations - 1;

    [x1, x2] = load_counts(mat, start, fin);

    [bin_v1, ~] = hist_with_errs(x1, edges);
    [bin_v2, ~] = hist_with_errs(x2, edges);

    [r, re] = ratio_and_err(bin_v1, bin_v2);

    sample_points_ratio(j-1, :)     = r;
    sample_points_ratio_err(j-1, :) = re;
end

%% Compute ratio of analytical to empirical uncertainty
sample_points_mean = mean(sample_points_ratio, 1);
sample_points_std  = std(sample_points_ratio, 0, 1);

moment_4 = moment(sample_points_ratio, 4);
sigma_b_2 = (moment_4 - (n_samples-2)/(n_samples) .* (var(sample_points_ratio,0,1).^2)) ...
            ./ n_samples ./ (4 .* var(sample_points_ratio,0,1));

% Alternative: Gaussian approximation for std error
% sigma_b   = sample_points_std ./ sqrt(2*(n_samples-1));
% sigma_b_2 = sigma_b.^2;

ratio_std     = ratio_bin_err ./ sample_points_std;
ratio_std_err = sqrt(sigma_a_2 ./ sample_points_std.^2 ...
              + sigma_b_2 .* ratio_bin_err.^2 ./ sample_points_std.^4);

print_ratios(ratio_std, ratio_std_err);

%% Plot: analytical vs empirical uncertainty
valid = isfinite(ratio_bin) & isfinite(ratio_bin_err) & ...
        isfinite(sample_points_mean) & isfinite(sample_points_std);

bin_c_p   = bin_c(valid);
ratio_p   = ratio_bin(valid);
ratio_e_p = ratio_bin_err(valid);
sp_mean_p = sample_points_mean(valid);
sp_std_p  = sample_points_std(valid);

figure_idx = figure_idx + 1;
figure(figure_idx); hold on; grid on; box on;

errorbar(bin_c_p, ratio_p,   ratio_e_p, '*');
errorbar(bin_c_p, sp_mean_p, sp_std_p,  '*', 'Color', 'red');
plot(bin_c, r_theory, '--', 'LineWidth', 1.2);
yline(1,            'LineStyle', '--');
yline(exp(epsilon),  '-.', 'Color', [1, 0.5, 0.3125]);
yline(exp(-epsilon), '-.', 'Color', [1, 0.5, 0.3125]);

legend('\sigma_{PLRV}', '\sigma_{sample}', 'Interpreter', 'tex');
title(["DP Definition - Privacy Loss r.v.", ...
       "Library: " + library + " | Bin: " + bin ...
       + " | $\epsilon$ = " + epsilon + " | N = " + iterations ...
       + " | Bin size = " + bin_size], ...
       'Interpreter', 'latex', 'FontSize', 12);
ylabel('Privacy Loss'); xlabel('Samples');
xlim([bin_c(1), bin_c(end)]);

% Overlay individual sample points
% for j = 1:size(sample_points_ratio, 1)
%     plot(bin_c, sample_points_ratio(j, :), '.', 'Color', 'green');
% end

%% Functions
function [x1, x2] = load_counts(mat, start, fin)
    x1 = mat.dp1(start:fin);
    x2 = mat.dp2(start:fin);
end

function [min_v, max_v, edges, centers] = make_edges_from_data(x1, x2, step, midpoint)
    min_v   = min([min(x1), min(x2)]);
    max_v   = max([max(x1), max(x2)]);
    c1      = midpoint - step/2;
    cs      = sort(c1:-step:min_v);
    centers = [cs, c1+step : step : max_v];
    edges   = (centers(1) - step/2) : step : (centers(end) + step/2);
end

function [bin_v, p_err] = hist_with_errs(x, edges)
    bin_v = histcounts(x, edges);
    [~, p_err] = Binomial_err(bin_v);
end

function [ratio_v, ratio_err] = ratio_and_err(v1, v2)
    n = numel(v1);
    ratio_v   = NaN(1, n);
    ratio_err = NaN(1, n);
    nz = (v1 ~= 0) & (v2 ~= 0);
    ratio_v(nz)   = v1(nz) ./ v2(nz);
    ratio_err(nz) = sqrt(v1(nz) ./ (v2(nz).^2) + (v1(nz).^2) ./ (v2(nz).^3));
end

function print_ratios(ratio_std, ratio_std_err)
    valid = isfinite(ratio_std) & isfinite(ratio_std_err) ...
          & (ratio_std ~= 0) & (ratio_std_err ~= 0);
    ratio_std     = ratio_std(valid);
    ratio_std_err = ratio_std_err(valid);

    n = numel(ratio_std);
    fprintf('%3d | % .6g ± %.6g\n', [(1:n).', ratio_std(:), ratio_std_err(:)]');

    k = 2;
    z = (ratio_std - 1) ./ ratio_std_err;
    count_2k = sum(z > k);

    p0   = 2 * (1 - normcdf(k));
    pval = binocdf(count_2k, n, p0, 'upper');

    fprintf('%d/%d = %.3f%% values above 2sigma (95%%), p-value = %.3g\n', ...
            count_2k, n, 100*count_2k/n, pval);
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
    if nargin < 8, xlim_override = [bin_c(1), bin_c(end)]; end
    figure_idx = figure_idx + 1;
    figure(figure_idx); hold on; grid on; box on;
    errorbar(bin_c, ratio_v, ratio_err, '*');
    plot(bin_c, r_theory, '--', 'LineWidth', 1.2);
    yline(1,            'LineStyle', '--');
    yline(exp(epsilon),  '-.', 'Color', [1, 0.5, 0.3125]);
    yline(exp(-epsilon), '-.', 'Color', [1, 0.5, 0.3125]);
    title(title_lines, 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Privacy Loss'); xlabel('Samples');
    xlim(xlim_override);
end