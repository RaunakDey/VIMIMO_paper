function out = welch_t_from_chains(chain1, chain2, n)
% WELCH_T_FROM_CHAINS  Welch two-sample t-test for two chains
% out = welch_t_from_chains(chain1, chain2)
% out = welch_t_from_chains(chain1, chain2, n)
%
% Inputs:
%  chain1, chain2 : vectors of samples (MCMC/posterior chains or any samples)
%  n              : optional number of draws to use from each chain
%                   (if omitted uses all available samples)
%                   if n > length(chain) sampling is with replacement
%
% Output (struct):
%  out.t         - Welch t statistic
%  out.df        - Welch-Satterthwaite degrees of freedom
%  out.p_two     - two-sided p-value
%  out.p_less    - one-sided p(mean1 < mean2)
%  out.p_greater - one-sided p(mean1 > mean2)
%  out.mean1, out.mean2
%  out.var1, out.var2
%  out.n1, out.n2 - sample sizes used

% validate and clean
if nargin < 3, useN = []; else useN = n; end
x = chain1(:); x = x(~isnan(x));
y = chain2(:); y = y(~isnan(y));
if isempty(x) || isempty(y)
    error('Both chains must contain at least one non-NaN sample.');
end

% draw samples if requested
if isempty(useN)
    xS = x;
    yS = y;
else
    if ~(isscalar(useN) && useN>0 && useN==floor(useN))
        error('n must be a positive integer scalar.');
    end
    nx = length(x); ny = length(y);
    if useN <= nx
        idx = randperm(nx, useN);
    else
        idx = randsample(nx, useN, true);
    end
    xS = x(idx);
    if useN <= ny
        idy = randperm(ny, useN);
    else
        idy = randsample(ny, useN, true);
    end
    yS = y(idy);
end

n1 = length(xS);
n2 = length(yS);
m1 = mean(xS);
m2 = mean(yS);
s1 = var(xS,0); % unbiased sample variance (div by n-1)
s2 = var(yS,0);

% standard error for difference in means (Welch)
se = sqrt(s1/n1 + s2/n2);

% t statistic
tstat = (m1 - m2) / se;

% Welch-Satterthwaite degrees of freedom
num = (s1/n1 + s2/n2)^2;
den = ( (s1^2) / ( (n1^2) * (n1-1) ) ) + ( (s2^2) / ( (n2^2) * (n2-1) ) );
df = num / den;

% guard
if ~isfinite(df) || df < 1
    df = 1;
end

% p-values
p_two = 2 * (1 - tcdf(abs(tstat), df));
p_less = tcdf(tstat, df);         % P(mean1 < mean2)
p_greater = 1 - tcdf(tstat, df); % P(mean1 > mean2)

% pack output
out.t = tstat;
out.df = df;
out.p_two = p_two;
out.p_less = p_less;
out.p_greater = p_greater;
out.mean1 = m1;
out.mean2 = m2;
out.var1 = s1;
out.var2 = s2;
out.n1 = n1;
out.n2 = n2;

res = out;

% Print summary
fprintf('Welch two-sample t-test:\n');
fprintf('  mean1 = %.4f, mean2 = %.4f\n', res.mean1, res.mean2);
fprintf('  t(%0.1f) = %.3f, two-sided p = %.4g\n', res.df, res.t, res.p_two);

% Decision rule (alpha = 0.05 by default)
alpha = 0.05;
if res.p_two < alpha
    fprintf('=> Difference detected (reject equal means at alpha=%.2f)\n', alpha);
else
    fprintf('=> No significant difference (cannot reject equality at alpha=%.2f)\n', alpha);
end

end