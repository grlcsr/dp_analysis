%% Parameters
clear; clc; close all;

figure_idx = 0;

libraries          = {'opendp', 'ibm', 'opendp_mod', 'ibm_original'};
bin_names          = {'rap2', 'openssl', 'xor', 'mgb_bad'};

iterations     = 1000000;

bin     = bin_names{1};
library = libraries{2};

switch lower(library)
    case 'opendp', tot_iterations = 100e6;
    case 'opendp_mod', tot_iterations = 100e6;
    case 'ibm',    tot_iterations = 91e6;
    case 'ibm_original', tot_iterations = 100e6;
    otherwise,     tot_iterations = 0;
end
tot_data = floor(tot_iterations/iterations);

epsilon        = 0.1;
bin_size       = 1;

dbsize1        = 10000;
dbsize2        = 10001;
midpoint = (dbsize1 + dbsize2)/2;

mat_name = "./data_mats/eps_" + sprintf('%.2f', epsilon) + "/" + library + "/" + bin + "/0H.mat";
%mat_name = "./data_mats/eps_" + sprintf('%.2f', epsilon) + "/" + library + "/geom.mat";

% --- Build nested output path: frozen_entropy/eps_xx.xx/<library>/<bin>/bin_size_<N>.mat
out_root   = "./statistical_significance/";
eps_folder = sprintf("eps_%0.2f", epsilon);
lib_folder = string(library);
bin_folder = string(bin);

out_dir     = fullfile(out_root, eps_folder, lib_folder, bin_folder);
out_mat_name = fullfile(out_dir, sprintf("bin_size_%d.mat", bin_size));

% Ensure folders exist
if ~isfolder(out_dir)
    mkdir(out_dir);
end

%% Compute metrics and save to .mat
mat      = load(mat_name);
[x1, x2] = load_counts(mat, 1, iterations);

% Common edges
[min_v, max_v, edges, bin_c] = make_edges_from_data(x1, x2, bin_size, midpoint);

% Histograms + binomial errors
[bin_v1, p_err1] = hist_with_errs(x1, edges);
[bin_v2, p_err2] = hist_with_errs(x2, edges);

% Ratio
[ratio_bin, ratio_bin_err] = ratio_and_err(bin_v1, bin_v2);

r_theory = exp(epsilon * sign(midpoint - bin_c));
pval = chi2_report(bin_c, ratio_bin, ratio_bin_err, r_theory);
%chi_report_68(bin_c, ratio_bin, ratio_bin_err, r_theory);

 [figure_idx] = plot_distributions(bin_c, bin_v1, p_err1, bin_v2, p_err2, ...
            ["DP Definition - Distributions", "No H manipulation - $\epsilon$ = " + epsilon + " - N = " + iterations + " - bin size = " + bin_size], ...
            figure_idx);

 [figure_idx] = plot_privacy_loss(bin_c, ratio_bin, ratio_bin_err, r_theory, ...
            ["DP Definition - Privacy Loss rv.", "No H manipulation - $\epsilon$ = " + epsilon + " - N = " + iterations + " - bin size = " + bin_size], ...
            figure_idx, epsilon, [bin_c(1), bin_c(end)]);
 
% Error propagation
sigma_a_2 = var(x1) * ((-4*bin_v1.*bin_v2 + 2*bin_v1*iterations + bin_v2*iterations) ./ (bin_v2.^3 * iterations)).^2 + ...
    var(x2) * (bin_v1 .* (4*bin_v1.*bin_v2 - 3*bin_v1*iterations - 2*bin_v2*iterations) ./ (bin_v2 .^4 * iterations)).^2;
%% Now read the sample points for significance
% Reading the sample points to check significance
sample_points_ratio = zeros(tot_data-1, length(ratio_bin));
sample_points_ratio_err = zeros(tot_data-1, length(ratio_bin));

for j = 2:tot_data
    start = 1 + (j - 1) * iterations;
    fin = start + iterations - 1;
    
    [x1, x2] = load_counts(mat, start, fin);
        
    % Histograms + binomial errors
    % Edges are the same as original histogram
    [bin_v1, p_err1] = hist_with_errs(x1, edges);
    [bin_v2, p_err2] = hist_with_errs(x2, edges);
    
    % Ratio + theory
    [ratio_bin_new, ratio_bin_err_new] = ratio_and_err(bin_v1, bin_v2);

    sample_points_ratio(j-1, :) = ratio_bin_new;
    sample_points_ratio_err(j-1, :) = ratio_bin_err_new;


    pval = chi2_report(bin_c, ratio_bin_new, ratio_bin_err_new, r_theory);
    %chi_report_68(bin_c, ratio_bin_new, ratio_bin_err_new, r_theory);

    if pval < 0.05
        fprintf('Significant result found for sample %d with p-value = %.4g\n', j, pval);
    end
end

save(out_mat_name, "bin_c", "edges", "ratio_bin", "ratio_bin_err", ...
     "sigma_a_2", "sample_points_ratio", "sample_points_ratio_err");

%% Load data and plot
load(out_mat_name, "bin_c", "edges", "ratio_bin", "ratio_bin_err", ...
     "sigma_a_2", "sample_points_ratio", "sample_points_ratio_err");

sample_points_mean = mean(sample_points_ratio, 1);
sample_points_std = std(sample_points_ratio, 0, 1);

moment_4 = moment(sample_points_ratio, 4);
sigma_b_2 = ((moment_4 - (tot_data-3)/(tot_data-1) .* (var(sample_points_ratio,0,1).^2)) ./ size(sample_points_ratio,1)) ./ (4 .* var(sample_points_ratio,0,1));

% Alternative for std error is gaussian approximation:
% n = tot_data - 1;
% sigma_b = sample_points_std ./ sqrt(2*(n-1));
% sigma_b_2 = sigma_b.^2;

ratio_std = ratio_bin_err./sample_points_std;
ratio_std_err = sqrt(sigma_a_2 ./ sample_points_std.^2 + sigma_b_2 .* ratio_bin_err.^2 ./ sample_points_std.^4);

print_ratios(ratio_std, ratio_std_err);

% Clean NaNs/Inf and align all series to the same valid indices
valid = isfinite(ratio_bin) & isfinite(ratio_bin_err) & ...
        isfinite(sample_points_mean) & isfinite(sample_points_std);

bin_c_p   = bin_c(valid);
ratio_p   = ratio_bin(valid);
ratio_e_p = ratio_bin_err(valid);
sp_mean_p = sample_points_mean(valid);
sp_std_p  = sample_points_std(valid);

% Plotting
midpoint = (dbsize1 + dbsize2)/2;
r_theory = exp(epsilon * sign(midpoint - bin_c));

figure_idx = figure_idx + 1;
figure(figure_idx); hold on; grid on; box on;

% errorbar(bin_c, ratio_bin, ratio_bin_err, '*');
% errorbar(bin_c, sample_points_mean, sample_points_std, '*', 'Color', 'red')
errorbar(bin_c_p, ratio_p,   ratio_e_p, '*');
errorbar(bin_c_p, sp_mean_p, sp_std_p,  '*', 'Color', 'red');
plot(bin_c, r_theory, '--', 'LineWidth', 1.2);
yline(1,  'LineStyle', '--'); 
yline(exp(epsilon),   '-.', 'Color', [1, 0.5, 0.3125]);
yline(exp(-epsilon),  '-.', 'Color', [1, 0.5, 0.3125]);
yline(1, 'LineStyle', '--');
legend('\sigma_{PLRV}', '\sigma_{sample}', 'Interpreter', 'tex')

title(["DP Definition - Privacy Loss r.v.", ...
       "Library: " + library + " | Bin: " + bin + " | $\epsilon$ = " + epsilon + " | N = "  + ...
       iterations + " | Bin size = " + bin_size], ...
       'Interpreter','latex','FontSize',12);

ylabel('Privacy Loss'); xlabel('Samples');
xlim([bin_c(1), bin_c(end)]);

% for j = 1:size(sample_points_ratio, 1)
%     plot(bin_c, sample_points_ratio(j, :), '.', 'Color', 'green');
% end

%% Sign test - both wrt r_theory and sample per sample
fprintf("Sign test analysis on samples from %s with binaries %s histbin size %d\n", library, bin, bin_size);
% Matlab function
idx_lower = bin_c_p < midpoint;
s = ones(size(bin_c_p)); s(~idx_lower) = -1;

dlog = log(ratio_p) - s .* epsilon;
valid = isfinite(dlog);
dlog = dlog(valid);
pval = signtest(dlog, 0);
fprintf('  Sign test (two tails, log-space): n=%d, p-value=%.4g\n', numel(dlog), pval);

% Option one: wrt r_theory not considering error bars
d = zeros(size(ratio_p));
d(idx_lower) = ratio_p(idx_lower) - exp(epsilon);
d(~idx_lower) = exp(-epsilon) - ratio_p(~idx_lower);

is_tie = ~isfinite(d) | (d == 0); % need to remove zeros (and infinites/NaNs in case)
signs  = sign(d(~is_tie));

k_pos  = sum(signs > 0);
n  = numel(signs);
k_min  = min(k_pos, n - k_pos);

pval  = min(2 * binocdf(k_min, n, 0.5), 1);
fprintf('  Sign test (two tails): n = %d, k_pos = %d, p-value = %.4g \n', n, k_pos, pval);

% Option 2: we consider error bars -> not really a sign test cuz we
% compare how many values are within the ksigma range with expected
% binomial
k = 3;
t = nan(size(ratio_p));
t(idx_lower)  = exp(epsilon);
t(~idx_lower) = exp(-epsilon);

z = (ratio_p - t) ./ ratio_e_p;

valid = isfinite(z) & isfinite(ratio_e_p) & (ratio_e_p > 0);
z = z(valid);

out = abs(z) > k;
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

s = ones(size(ratio_p)); s(~idx_lower) = -1;
dlog = log(ratio_p) - s.*epsilon;

% delta-method SE for log(r): sigma_log ≈ sigma_ratio / ratio
sigma_log = ratio_e_p ./ ratio_p;

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


%% Functions
function [x1, x2] = load_counts(mat, start, fin)
    x1 = mat.dp1(start:fin);
    x2 = mat.dp2(start:fin);
end

function [min_v, max_v, edges, centers] = make_edges_from_data(x1, x2, step, midpoint)
    min_v   = min(min(x1, x2));
    max_v   = max(max(x1, x2));
    c1      = midpoint - step/2;
    cs      = sort(c1:-step:min_v);
    centers = [cs, c1 + step:step:max_v];
    edges   = centers(1) - step/2 : step : centers(end) + step/2;
    % centers = min_v:step:max_v;
    % edges   = (min_v - 0.5*step) : step : (max_v + 0.5*step);
end

function [bin_v, p_err] = hist_with_errs(x, edges)
    bin_v = histcounts(x, edges);
    [~, p_err] = Binomial_err(bin_v);
end

function [ratio_v, ratio_err] = ratio_and_err(v1, v2)
    n = numel(v1);
    ratio_v  = NaN(1, n);
    ratio_err = NaN(1, n);
    nz = (v1 ~= 0) & (v2 ~= 0);
    ratio_v(nz) = v1(nz) ./ v2(nz);
    ratio_err(nz) = sqrt( v1(nz) ./ (v2(nz).^2) + (v1(nz).^2) ./ (v2(nz).^3) );
end

function print_ratios(ratio_std, ratio_std_err)
    valid = (ratio_std ~= 0) & (ratio_std_err ~= 0) & isfinite(ratio_std) & isfinite(ratio_std_err);
    ratio_std = ratio_std(valid); 
    ratio_std_err = ratio_std_err(valid);

    n = numel(ratio_std);
    fprintf('%3d | % .6g ± %.6g\n', [ (1:n).', ratio_std(:), ratio_std_err(:) ]');
    
    k = 2;
    z = (ratio_std - 1)./ratio_std_err;
    count_2k = length(z(z>k));

    p0 = 2 * (1 - normcdf(k));
    pval = binocdf(count_2k, n, p0, 'upper');

    fprintf('%d/%d = %3f%% values above 2sigma (95%%), p-value = %.3g \n', count_2k, length(z), 100*count_2k/length(z), pval);

end

function pval = chi2_report(bin_c, ratio_bin, ratio_bin_err, r_theory)
    ndf = 0;
    chi2_fit = 0;
    
    for j = 1:length(bin_c)
        if(ratio_bin_err(j) ~= 0 && ~isnan(ratio_bin_err(j)) && ~isnan(ratio_bin(j)) && ratio_bin(j) ~= 0)
            chi2_fit = chi2_fit + ((r_theory(j)-ratio_bin(j))/ratio_bin_err(j))^2;
            ndf = ndf +1;
        end
    end

    % normalize
    chi2_red = chi2_fit/ndf;
    chi2_red_std = sqrt(2/ndf);

    pval = 1 - chi2cdf(chi2_fit, ndf-1);

    fprintf("  χ² fit with theoretical epsilon:\n");
    fprintf("  χ²_red = %.5f ± %.5f   (ndf = %d), p-value = %.3g\n", chi2_red, chi2_red_std, ndf, pval);
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