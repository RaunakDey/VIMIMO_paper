function out = ttest_from_chains(chain1, chain2, n)
% TTEST_FROM_CHAINS  Welch two-sample t-test from MCMC/posterior chains
% out = ttest_from_chains(chain1, chain2)
% out = ttest_from_chains(chain1, chain2, n)
%
% Inputs
%  chain1, chain2 : vectors of samples (column or row)
%  n              : optional number of observations to draw from each chain
%                   if omitted uses full lengths. If n > length(chain),
%                   sampling is with replacement.
%
% Output (struct)
%  out.t        : t statistic (Welch)
%  out.df       : Welch-Satterthwaite degrees of freedom
%  out.p_two    : two-sided p-value
%  out.p_less   : P(mean1 < mean2) using t-distribution (one-sided)
%  out.p_greater: P(mean1 > mean2) using t-distribution (one-sided)
%  out.meanDiff : mean(chain1_sample) - mean(chain2_sample)
%  out.n1, out.n2: sample sizes used
%
% Notes
%  - Handles NaNs by removing them before sampling.
%  - Uses sampling without replacement when possible, otherwise with replacement.

if nargin < 3
    useAll = true;
else
    useAll = false;
    if ~(isscalar(n) && n>0 && n==floor(n))
        error('n must be a positive integer scalar.');
    end
end

% enforce vectors and remove NaNs
chain1 = chain1(:);
chain2 = chain2(:);
chain1 = chain1(~isnan(chain1));
chain2 = chain2(~isnan(chain2));

if isempty(chain1) || isempty(chain2)
    error('Both chains must contain at least one non-NaN sample.');
end

% sample n if requested
if useAll
    x = chain1;
    y = chain2;
else
    n1 = length(chain1);
    n2 = length(chain2);
    if n <= n1
        idx1 = randperm(n1,n);
    else
        idx1 = randsample(n1,n,true);
    end
    if n <= n2
        idx2 = randperm(n2,n);
    else
        idx2 = randsample(n2,n,true);
    end
    x = chain1(idx1);
    y = chain2(idx2);
end

n1 = length(x);
n2 = length(y);
mx = mean(x);
my = mean(y);
sx2 = var(x,1); % population var (divide by n) consistent with formula below
sy2 = var(y,1);

% Use sample variances (unbiased) in DF formula requires sample var with (n-1)
% so compute unbiased sample variances:
sx2_unb = var(x,0); % divides by (n-1)
sy2_unb = var(y,0);

% t-statistic (Welch)
se = sqrt(sx2./n1 + sy2./n2); % note sx2 used is population var; equivalent to unbiased for t formula
tstat = (mx - my) / se;

% Welch-Satterthwaite degrees of freedom
% using population variances in numerator but unbiased in denom; equivalent standard formula:
num = (sx2./n1 + sy2./n2).^2;
den = ( (sx2.^2) ./ ( (n1.^2) .* (n1-1) ) ) + ( (sy2.^2) ./ ( (n2.^2) .* (n2-1) ) );
df = num ./ den;

% two-sided p-value
p_two = 2 * (1 - tcdf(abs(tstat), df));

% one-sided p-values (using t-dist)
p_less = tcdf(tstat, df);          % P(T <= t) => evidence mean1 < mean2
p_greater = 1 - tcdf(tstat, df);  % P(T >= t) => evidence mean1 > mean2

% Package output
out.t = tstat;
out.df = df;
out.p_two = p_two;
out.p_less = p_less;
out.p_greater = p_greater;
out.meanDiff = mx - my;
out.n1 = n1;
out.n2 = n2;

end