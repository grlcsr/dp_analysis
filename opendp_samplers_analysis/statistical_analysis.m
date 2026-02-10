clc; close all;

lap_scale_name = "10";          % Laplace scale (b = 1/epsilon), here b is an integer
lap_scale      = 10;
matfile        = "samples_" + lap_scale_name + "_bad.mat";

% I/O: load once and pass the struct around 
S = load(matfile);

% -------------------------------------------------------------------------
% Uniformity tests (don't depend on lap_scale)
% -------------------------------------------------------------------------
test_uniformity(S);             % unif last digit + plot
test_uniformity_unif_below(S);  % unif_below uniformity + plot

% -------------------------------------------------------------------------
% Bernoulli tests
% -------------------------------------------------------------------------
test_bernoulli_standard(S);            % tests S.bernoulli vs p = 0.5
test_bernoulli_rat(S, lap_scale);      % tests S.bernoulli_rat vs p = 1/lap_scale
test_bernoulli_exp1(S, lap_scale);     % tests over U: E[exp(-U/b)]
test_bernoulli_exp(S, lap_scale);      % sme as bafore

% -------------------------------------------------------------------------
% Geometric (slow vs fast CKS20)
% -------------------------------------------------------------------------
test_geom_slow(S);                     % x = 1
test_geom_fast(S, lap_scale);          % x = 1/lap_scale

% -------------------------------------------------------------------------
% Discrete Laplace
% -------------------------------------------------------------------------
test_dlap(S, lap_scale);

%% ========================================================================
%% Uniformity tests
%% ========================================================================

function test_uniformity(S)
    % Quick sanity check on S.unif: last decimal digit should be ~uniform on 0..9

    u = S.unif;

    if ~isstring(u) && ~iscellstr(u)
        error('unif is not a string array or cellstr; check import.');
    end
    if iscell(u)
        u = string(u);
    end

    n_u = numel(u);
    lastDigit = zeros(n_u, 1);

    for i = 1:n_u
        s = char(u(i));
        c = s(end);  % last character
        d = double(c) - double('0');
        if d < 0 || d > 9
            error('Found non-digit last char in unif at row %d: "%s"', i, s);
        end
        lastDigit(i) = d;
    end

    edges   = -0.5:9.5;
    counts  = histcounts(lastDigit, edges);
    expected_u = n_u / 10;

    chi2_u = sum((counts - expected_u).^2 ./ expected_u);
    p_u    = 1 - chi2cdf(chi2_u, 9);

    fprintf('=== unif last-decimal-digit uniformity (0-9) ===\n');
    fprintf('n       = %d\n', n_u);
    fprintf('counts  ='); fprintf(' %d', counts); fprintf('\n');
    fprintf('chi2    = %.4g (df = 9)\n', chi2_u);
    fprintf('p-value = %.4g\n\n', p_u);

    cats  = 0:9;
    p0    = 1/10;
    sigma = sqrt(n_u * p0 * (1 - p0));  % rough binomial sigma, same for all bins

    figure;
    hold on;
    errorbar(cats, counts, sigma * ones(size(cats)), '*', 'LineWidth', 1.2);
    plot(cats, expected_u * ones(size(cats)), '--', 'LineWidth', 1.2);
    hold off;
    xlabel('Last decimal digit of unif');
    ylabel('Counts');
    title(sprintf('Uniformity last-digit: \\chi^2 = %.3g (dof = 9), p = %.3g', chi2_u, p_u), 'Interpreter', 'tex', 'FontSize', 15);
    legend('Observed counts', 'Expected', 'Location', 'best');
    grid on;
end

function test_uniformity_unif_below(S)
    % Chi-square uniformity test for the integer-valued unif_below

    x = S.unif_below;

    if ~isnumeric(x)
        error('unif_below is not numeric; check import in make_mat_from_log.');
    end

    x = x(:);
    x = x(~isnan(x));  % drop NaNs from missing rows

    if isempty(x)
        warning('test_uniformity_unif_below: no finite samples left, skipping.');
        return;
    end

    [support, ~, idx] = unique(x);
    M = numel(support);
    counts = accumarray(idx, 1, [M 1]);

    n = numel(x);
    expected = n / M;

    chi2 = sum((counts - expected).^2 ./ expected);
    df   = M - 1;
    p    = 1 - chi2cdf(chi2, df);

    fprintf('=== unif_below chi-square uniformity over %d categories ===\n', M);
    fprintf('support ='); fprintf(' %g', support); fprintf('\n');
    fprintf('counts  ='); fprintf(' %d', counts); fprintf('\n');
    fprintf('chi2    = %.4g (df = %d)\n', chi2, df);
    fprintf('p-value = %.4g\n\n', p);

    p_cat = 1 / M;
    sigma = sqrt(n * p_cat * (1 - p_cat));  % again same-ish sigma for all

    figure;
    hold on;
    errorbar(support, counts, sigma * ones(size(support)), '*', 'LineWidth', 1.2);
    plot(support, expected * ones(size(support)), '--', 'LineWidth', 1.2);
    hold off;
    xlabel('unif\_below value');
    ylabel('Counts');
    title(sprintf('unif\\_below: \\chi^2 = %.3g (dof = %d), p = %.3g', chi2, df, p), 'Interpreter', 'tex', 'FontSize', 15);
    legend('Observed counts', 'Expected', 'Location', 'best');
    grid on;
end

%% ========================================================================
%% Bernoulli tests
%% ========================================================================

function test_bernoulli_generic(x, p_theory, name)
    % Generic Bernoulli test: compares empirical p_hat with p_theory.
    % Nothing fancy: z-test + chi^2 on 0/1 counts.

    x = x(:);

    if isnumeric(x)
        x = x(~isnan(x));  % drop NaNs from missing CSV cells
    end

    if isempty(x)
        warning('test_bernoulli_generic(%s): no finite samples, skipping.\n', name);
        return;
    end

    if ~islogical(x)
        x = x ~= 0;  % coerce 0/1 -> logical
    end

    n  = numel(x);
    k1 = sum(x);
    k0 = n - k1;

    p_hat = k1 / n;

    se = sqrt(p_theory * (1 - p_theory) / n);
    z  = (p_hat - p_theory) / se;
    p_z = 2 * normcdf(-abs(z), 0, 1);

    expected1 = n * p_theory;
    expected0 = n * (1 - p_theory);
    chi2 = (k1 - expected1)^2 / expected1 + (k0 - expected0)^2 / expected0;
    df   = 1;
    p_chi = 1 - chi2cdf(chi2, df);

    fprintf('=== %s Bernoulli test ===\n', name);
    fprintf('n        = %d\n', n);
    fprintf('k1       = %d (p_hat = %.6f)\n', k1, p_hat);
    fprintf('p_theory = %.6f\n', p_theory);
    fprintf('z-test:  z = %.4f,  p = %.4g\n', z, p_z);
    fprintf('chi^2:   chi2 = %.4f (df = %d), p = %.4g\n\n', chi2, df, p_chi);

    cats     = [0 1];
    counts   = [k0 k1];
    expected = [expected0 expected1];

    sigma = sqrt(n * p_theory * (1 - p_theory));

    figure;
    hold on;
    errorbar(cats, counts, sigma * ones(size(cats)), '*', 'LineWidth', 1.2);
    plot(cats, expected, 'x', 'LineWidth', 1.2, 'MarkerSize', 12);
    hold off;
    xlim([-0.5 1.5]);
    xticks([0 1]);
    xlabel(sprintf('%s value', name));
    ylabel('Counts');
    title(sprintf('%s: p_{hat}=%.4f, p_{th}=%.4f, \\chi^2=%.3g, p=%.3g', ...
                  name, p_hat, p_theory, chi2, p_chi), 'FontSize', 15, 'Interpreter','tex');
    legend('Observed counts', 'Expected', 'Location', 'best');
    grid on;
end

function test_bernoulli_standard(S)
    % S.bernoulli vs p = 0.5

    x = S.bernoulli;
    p_theory = 0.5;
    test_bernoulli_generic(x, p_theory, 'bernoulli');
end

function test_bernoulli_rat(S, lap_scale)
    % bernoulli_rat should be Bernoulli(1/lap_scale) (i.e., epsilon)

    if nargin < 2 || isempty(lap_scale)
        error('test_bernoulli_rat: lap_scale must be provided (b = 1/epsilon).');
    end

    x = S.bernoulli_rat;

    p_theory = 1 / lap_scale;   % epsilon
    test_bernoulli_generic(x, p_theory, 'bernoulli\_rat');
end

function p_mix = bern_exp_mixture_p(lap_scale)
    % E[exp(-U/b)] with U ~ Unif{0..b-1}, b = lap_scale
    % Just the closed-form sum of a finite geometric series.

    b = lap_scale;
    r = exp(-1 / b);
    p_mix = (1 / b) * (1 - exp(-1)) / (1 - r);
end

function test_bernoulli_exp1(S, lap_scale)
    % Test S.bernoulli_exp1 against the mixture p = E[exp(-U/b)]

    if nargin < 2 || isempty(lap_scale)
        error('test_bernoulli_exp1: lap_scale must be provided (b = 1/epsilon).');
    end

    x = S.bernoulli_exp1;
    p_theory = bern_exp_mixture_p(lap_scale);

    test_bernoulli_generic(x, p_theory, 'bernoulli\_exp1');
end

function test_bernoulli_exp(S, lap_scale)
    % Same mixture as above, but for S.bernoulli_exp

    if nargin < 2 || isempty(lap_scale)
        error('test_bernoulli_exp: lap_scale must be provided.');
    end

    x = S.bernoulli_exp;
    p_theory = bern_exp_mixture_p(lap_scale);

    test_bernoulli_generic(x, p_theory, 'bernoulli\_exp');
end

%% ========================================================================
%% Geometric tests
%% ========================================================================

function test_geometric_generic(x, x_param, name)
    % Check geometric samples vs theory for given x_param.
    % P(K = k) = (1 - e^{-x_param}) e^{-x_param k},  k = 0,1,...

    x = double(x(:));
    x = x(~isnan(x));     % drop NaNs from missing CSV cells

    if isempty(x)
        warning('test_geometric_generic(%s): no finite samples, skipping.\n', name);
        return;
    end

    if any(x < 0)
        error('%s: geometric samples must be >= 0', name);
    end

    n = numel(x);

    q = exp(-x_param);
    p = 1 - q;

    mean_th = 1 / (exp(x_param) - 1);   % = q/p
    var_th  = q / p^2;

    mean_hat = mean(x);
    var_hat  = var(x, 1);   % population variance

    se_mean = sqrt(var_th / n);
    z_mean  = (mean_hat - mean_th) / se_mean;
    p_mean  = 2 * normcdf(-abs(z_mean), 0, 1);

    fprintf('=== %s geometric test ===\n', name);
    fprintf('n         = %d\n', n);
    fprintf('mean_hat  = %.6f,  mean_th = %.6f\n', mean_hat, mean_th);
    fprintf('var_hat   = %.6f,  var_th  = %.6f\n', var_hat,  var_th);
    fprintf('mean z    = %.4f,  p = %.4g\n\n', z_mean, p_mean);

    % Histogram vs theoretical pmf (truncated to something not insane)
    K_max_data = max(x);
    K_plot = K_max_data; %min(K_max_data, 100);  % arbitrary but fine

    edges  = -0.5:1:(K_plot + 0.5);
    counts = histcounts(x, edges);
    ks = 0:K_plot;

    pmf_full   = p * q.^ks;
    mass_trunc = sum(pmf_full);
    pmf_trunc  = pmf_full / mass_trunc;
    expected   = n * pmf_trunc;

    % Binomial errors on the observed counts
    [~, sigma] = Binomial_err(counts);

    chi2 = sum((counts - expected).^2 ./ max(expected, eps));
    df   = numel(ks) - 1;
    p_chi = 1 - chi2cdf(chi2, df);

    fprintf('%s chi-square over k = 0..%d:\n', name, K_plot);
    fprintf('chi2 = %.4g (df = %d), p = %.4g\n\n', chi2, df, p_chi);

    figure;
    hold on;
    errorbar(ks, counts, sigma, '*', 'LineStyle', 'none', 'LineWidth', 1.2);
    plot(ks, expected, 'LineWidth', 1.5);
    hold off;

    xlabel(sprintf('%s value k', name));
    ylabel('Counts');
    title(sprintf('%s: mean z=%.2f (p=%.2g), \\chi^2=%.2g (p=%.2g)', ...
          name, z_mean, p_mean, chi2, p_chi), 'FontSize', 15, 'Interpreter','tex');
    legend('Observed', 'Expected', 'Location', 'best');
    grid on;
end

function test_geom_slow(S)
    % S.geom_slow, which should correspond to x_param = 1

    x = S.geom_slow;
    x_param = 1;
    test_geometric_generic(x, x_param, 'geom\_slow');
end

function test_geom_fast(S, lap_scale)
    % S.geom_fast, CKS20-style geometric with x_param = 1/lap_scale

    if nargin < 2 || isempty(lap_scale)
        error('test_geom_fast: lap_scale must be provided.');
    end

    x = S.geom_fast;
    x_param = 1 / lap_scale;
    test_geometric_generic(x, x_param, 'geom\_fast');
end

%% ========================================================================
%% Discrete Laplace tests
%% ========================================================================

function test_dlap(S, lap_scale)
    % Check S.dlap against discrete Laplace with scale = lap_scale

    if nargin < 2 || isempty(lap_scale)
        error('test_dlap: lap_scale must be provided.');
    end

    x = double(S.dlap(:));
    x = x(~isnan(x));      % drop NaNs from missing CSV cells

    if isempty(x)
        warning('test_dlap: no finite samples, skipping.');
        return;
    end

    n = numel(x);

    eps_val = 1 / lap_scale;
    q = exp(-eps_val);
    c = (1 - q) / (1 + q);

    mean_th = 0;
    var_th  = 2 * q / (1 - q)^2;

    mean_hat = mean(x);
    var_hat  = var(x, 1);

    se_mean = sqrt(var_th / n);
    z_mean  = (mean_hat - mean_th) / se_mean;
    p_mean  = 2 * normcdf(-abs(z_mean), 0, 1);

    fprintf('=== discrete Laplace test ===\n');
    fprintf('n         = %d\n', n);
    fprintf('mean_hat  = %.6f,  mean_th = %.6f\n', mean_hat, mean_th);
    fprintf('var_hat   = %.6f,  var_th  = %.6f\n', var_hat,  var_th);
    fprintf('mean z    = %.4f,  p = %.4g\n\n', z_mean, p_mean);

    % Now the histogram / chi-square stuff
    abs_x = abs(x);
    K_max_data = max(abs_x);
    K = K_max_data;

    support = -K:K;
    edges   = (support - 0.5);
    edges   = [edges, support(end) + 0.5];

    counts = histcounts(x, edges);

    pmf_full   = c * q.^abs(support);
    mass_trunc = sum(pmf_full);
    pmf_trunc  = pmf_full / mass_trunc;
    expected   = n * pmf_trunc;

    % Binomial errors for each count bin
    [~, sigma] = Binomial_err(counts);

    chi2 = sum((counts - expected).^2 ./ max(expected, eps));
    df   = numel(support) - 1;
    p_chi = 1 - chi2cdf(chi2, df);

    fprintf('discrete Laplace chi-square over x = -%d..%d:\n', K, K);
    fprintf('chi2 = %.4g (df = %d), p = %.4g\n\n', chi2, df, p_chi);

    % I like to look at normalized counts here 
    figure;
    hold on;
    errorbar(support, counts / n, sigma / n, '*');
    plot(support, expected / n, 'LineWidth', 1.5);
    hold off;

    xlabel('dlap value x');
    ylabel('Normalized counts');
    title(sprintf('discrete lap: mean z=%.2f (p=%.2g), \\chi^2=%.2g (p=%.2g)', ...
          z_mean, p_mean, chi2, p_chi), 'FontSize', 15, 'Interpreter','tex');
    legend('Observed', 'Expected', 'Location', 'best');
    grid on;
end


function [weight, p_err] = Binomial_err(x)
    % Rough binomial error estimate for histogram counts.
    % x      : vector of counts per bin
    % weight : 1 / sigma^2 per bin (0 if count is 0)
    % p_err  : sigma in counts for each bin

    x   = x(:)';      % row vector
    n   = length(x);
    Tot = sum(x);

    weight = zeros(1, n);
    p_err  = zeros(1, n);

    if Tot == 0
        return;
    end

    for j = 1:n
        if x(j) == 0
            weight(j) = 0;
            p_err(j)  = 0;
        else
            p        = x(j) / Tot;
            p_err(j) = sqrt(Tot * p * (1 - p));   % sigma in counts
            weight(j)= 1 / p_err(j)^2;
        end
    end
end
