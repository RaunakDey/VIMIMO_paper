function out = bayes_eq_test(chainA, chainB, varargin)
% BAYES_EQ_TEST  Bayesian equality test via ROPE probability
% out = bayes_eq_test(chainA, chainB)
% out = bayes_eq_test(chainA, chainB, 'rope', 0.1, 'nsamp', 10000, 'mode','abs', 'alpha_eq', 0.95)
%
% Inputs (name-value):
%  'rope'    : ROPE half-width (default: 0.1 * pooled std if mode='abs',
%              or 0.1 when mode='std' means 0.1 SDs). 
%  'mode'    : 'abs' (absolute ROPE) or 'std' (ROPE in units of pooled SD). Default 'abs'.
%  'nsamp'   : number of paired draws to form delta (default: use max length available; resample with replacement).
%  'cred'    : credible level for CrI, default 0.95.
%  'alpha_eq': threshold for declaring "same" (default 0.95). 'alpha_diff' uses 1-alpha_eq for "different".
%  'print'   : true/false (default true) to print human readable result.
%
% Output:
%  out.meanDiff, out.CrI, out.P_greater, out.P_inROPE (p_eq), out.ROPE, out.nsamp, out.decision

% Parse inputs
p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'rope',[], @(x) isnumeric(x) && isscalar(x));
addParameter(p,'mode','abs', @(s) any(strcmp(s,{'abs','std'})));
addParameter(p,'nsamp',[], @(x) isempty(x) || (isscalar(x) && x>0 && x==floor(x)));
addParameter(p,'cred',0.95, @(x) isnumeric(x) && x>0 && x<1);
addParameter(p,'alpha_eq',0.95, @(x) isnumeric(x) && x>0 && x<1);
addParameter(p,'print',true, @islogical);
parse(p,varargin{:});
opts = p.Results;

% Clean inputs
A = chainA(:); A = A(~isnan(A));
B = chainB(:); B = B(~isnan(B));
if isempty(A) || isempty(B)
    error('Both chains must contain at least one non-NaN sample.');
end

% Determine nsamp (use max available to reduce sampling variability)
if isempty(opts.nsamp)
    nsamp = max(length(A), length(B));
else
    nsamp = opts.nsamp;
end

% Draw paired samples (resample with replacement if needed)
if nsamp <= length(A)
    idxA = randperm(length(A), nsamp);
else
    idxA = randsample(length(A), nsamp, true);
end
if nsamp <= length(B)
    idxB = randperm(length(B), nsamp);
else
    idxB = randsample(length(B), nsamp, true);
end
As = A(idxA);
Bs = B(idxB);

% Compute delta
delta = As - Bs;

% Determine ROPE
pooledSD = sqrt( var(A,0)* (length(A)-1) + var(B,0)*(length(B)-1) ) ...
           / sqrt( (length(A)+length(B)-2) ); % pooled sample SD (unbiased)
if strcmp(opts.mode,'std')
    if isempty(opts.rope)
        rope = 0.1 * pooledSD; % default 0.1 SD
    else
        rope = opts.rope * pooledSD; % user provided factor in SD units
    end
else % mode 'abs'
    if isempty(opts.rope)
        % default absolute ROPE: 0.05 * range of pooled samples or 0.1 * pooledSD if range too small
        pooledRange = max([A;B]) - min([A;B]);
        if pooledRange > 0
            rope = 0.05 * pooledRange;
        else
            rope = 0.1 * pooledSD;
        end
    else
        rope = opts.rope;
    end
end

% Summaries
meanDiff = mean(delta);
lower = quantile(delta, (1-opts.cred)/2);
upper = quantile(delta, 1-(1-opts.cred)/2);
P_greater = mean(delta > 0);          % P(A > B)
P_inROPE = mean(abs(delta) < rope);   % Bayesian "p_eq"

% Decision rule
alpha_eq = opts.alpha_eq;
alpha_diff = 1 - alpha_eq;
if P_inROPE >= alpha_eq
    decision = 'same';
elseif P_inROPE <= alpha_diff
    decision = 'different';
else
    decision = 'inconclusive';
end

% Pack output
out.meanDiff = meanDiff;
out.CrI = [lower, upper];
out.P_greater = P_greater;
out.P_inROPE = P_inROPE;
out.ROPE = rope;
out.nsamp = nsamp;
out.pooledSD = pooledSD;
out.decision = decision;

% Print
if opts.print
    fprintf('Bayesian equality test (ROPE method)\n');
    fprintf(' nsamp = %d, ROPE = %.6g\n', nsamp, rope);
    fprintf(' meanDiff = %.6g, %.1f%% CrI = [%.6g, %.6g]\n', meanDiff, 100*opts.cred, lower, upper);
    fprintf(' P(A > B) = %.4f\n', P_greater);
    fprintf(' p_eq = P(|Δ| < ROPE) = %.4f\n', P_inROPE);
    if strcmp(decision,'same')
        fprintf(' => DECISION: SAME (p_eq >= %.3f)\n', alpha_eq);
    elseif strcmp(decision,'different')
        fprintf(' => DECISION: DIFFERENT (p_eq <= %.3f)\n', alpha_diff);
    else
        fprintf(' => DECISION: INCONCLUSIVE (%.3f < p_eq < %.3f)\n', alpha_diff, alpha_eq);
    end
    fprintf(' Note: this is Bayesian; p_eq is posterior probability of practical equivalence.\n\n');
end
end