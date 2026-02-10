%% Parameters
clear; clc; close all;

entropies  = {"M4", "M8", "M16", "M32", "M64"};

bits = {"clear", "set", "flip", "alternate", "correlation"};

libraries = {'opendp_mod','opendp', 'ibm'};
bin_names = {'rap2', 'openssl', 'xor', 'mgb_bad'};

iterations     = 1000000;

bin            = bin_names{1};
library        = libraries{5};

epsilon        = 0.1;
bin_sizes      = 5;
hs             = 1:length(entropies);
bs             = 1:length(bits);

dbsize1        = 10000;
dbsize2        = 10001;
midpoint = (dbsize1 + dbsize2)/2;

mat_name = "./data_mats/eps_" + sprintf('%.2f', epsilon) + "/" + library + "/0H.mat";
mat     = load(mat_name);

x1 = mat.dp1 - dbsize1;
x2 = mat.dp2 - dbsize2;
conc_x1 = cat(1, x1, x2);
step = 1;

pvals = []; 
low_pval = 0;
for j = 1:100
    start = 1 + (j - 1) * iterations;
    fin = start + iterations - 1;
    x1 = conc_x1(start:fin);

    min_v = min(x1); max_v = max(x1);
    bin_c = min_v:step:max_v;
    edges   = (min_v - 0.5*step) : step : (max_v + 0.5*step);
    [bin_v, bin_err] = hist_with_errs(x1, edges);
    
    bin_v_norm = bin_v/sum(bin_v);
    pError_norm = bin_err/sum(bin_v);
    scale = 0.1;
    
    vals = discrete_laplace_pdf(scale, bin_c);
        
    pval = chi2_report(bin_c, bin_v_norm, pError_norm, vals);

    if pval < 0.05
        low_pval = low_pval + 1;
    end
    pvals = [pvals, pval];
end

cdf_pvals = sort(pvals);
idx = (1:length(cdf_pvals)) * 1/length(pvals);

%% Fit

[xData, yData] = prepareCurveData(idx, cdf_pvals);

ft = fittype('poly1');
opts = fitoptions('Method', 'LinearLeastSquares');
opts.Robust = 'Bisquare';

% Fit model to data
[fitresult, gof] = fit(xData, yData, ft, opts);
conf = confint(fitresult, 0.68);

figure();
h = plot(fitresult, xData, yData , '*');
legend(h, 'Sample CDF', 'Fit', 'Location', 'Best', 'Interpreter', 'none');

xlabel('index', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('CDF', 'Interpreter', 'latex', 'FontSize', 13);
ylim([0, 1]);
grid on
title('Sample p-values CDF', 'Interpreter', 'latex', 'FontSize', 15);
info = strcat("a = ", num2str(fitresult.p1), " $\pm$ ", num2str(abs(conf(1,1) - fitresult.p1)));
text(0.1, 0.9, info, 'FontSize', 11, 'Interpreter', 'latex');
info = strcat("b = ", num2str(fitresult.p2), " $\pm$ ", num2str(abs(conf(2,2) - fitresult.p2)));
text(0.1, 0.85, info, 'FontSize', 11, 'Interpreter', 'latex');

[h,p,ksstat,cv] = kstest(pvals, 'CDF', makedist('Uniform',0,1));
fprintf('h = %d, p = %.4g, ksstat = %.4g, cv ≈ %.4g\n', h, p, ksstat, cv);


%% Functions 

function vals = discrete_laplace_pdf(scale, x)
    p = exp(-scale);
    vals = (1-p)/(1+p) * p.^(abs(x));
end

function [bin_v, p_err] = hist_with_errs(x, edges)
    bin_v = histcounts(x, edges);
    [~, p_err] = Binomial_err(bin_v);
end

function pval = chi2_report(bin_c, ratio_bin, ratio_bin_err, r_theory)
    ndf = 0;
    chi2_fit = 0;
    
    for j = 1:length(bin_c)
        if(ratio_bin_err(j) ~= 0)
            chi2_fit = chi2_fit + ((r_theory(j)-ratio_bin(j))/ratio_bin_err(j))^2;
            ndf = ndf +1;
        end
    end

    % normalize
    chi2_red = chi2_fit/ndf;
    chi2_red_std = sqrt(2/ndf);

    pval = 1 - chi2cdf(chi2_fit, ndf);

    fprintf("χ² fit with theoretical epsilon:\n");
    fprintf("  χ²_red = %.5f ± %.5f   (ndf = %d)\n", chi2_red, chi2_red_std, ndf);
    fprintf("  p-value = %.3g\n\n", pval);
end
