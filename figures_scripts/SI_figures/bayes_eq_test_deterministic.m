function out = bayes_eq_test_deterministic(chainA, chainB, varargin)
% BAYES_EQ_TEST_DETERMINISTIC  Bayesian ROPE equality test, deterministic sampling
% out = bayes_eq_test_deterministic(chainA, chainB)
% out = bayes_eq_test_deterministic(chainA, chainB, 'nsamp', 20000, ...
%       'mode','std','rope',0.1,'sampling','deterministic','cred',0.95)
%
% Behavior (no randomness by default):
%  - If length(chainA)==length(chainB) and nsamp is empty -> paired elementwise Δ = A - B
%  - Otherwise indices are chosen deterministically via linspace to get nsamp samples
%  - 'sampling' can be: 'paired','truncate','deterministic','resample'
%      'paired'    : requires equal lengths, uses A-B elementwise
%      'truncate'  : take first min(Na,Nb) pairs
%      'deterministic' : pick nsamp indices via linspace (default for unequal lengths)
%      'resample'  : random resampling with replacement (only if user requests)
%
% Options (name-value):
%  'nsamp'   : number of paired draws to form delta (default: [] means use full if equal,
%              else min(lengths)).
%  'mode'    : 'abs' or 'std' for ROPE specification (default 'abs')
%  'rope'    : ROPE half-width (scalar). If mode='std' rope is fraction of pooled SD.
%  'sampling': see above (default '')
%  'cred'    : credible level for CrI (default 0.95)
%  'alpha_eq': threshold for declaring "same" (default 0.95)
%  'print'   : true/false (default true)
%
% Output:
%  out.meanDiff, out.CrI, out.P_greater, out.P_inROPE, out.ROPE, out.nsamp, out.decision

% Parse inputs
p = inputParser;
addParameter(p,'nsamp',[], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0 && x==floor(x)));
addParameter(p,'mode','abs', @(s) any(strcmp(s,{'abs','std'})));
addParameter(p,'rope',[], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p,'sampling','', @(s) ischar(s) || isstring(s));
addParameter(p,'cred',0.95, @(x) isnumeric(x) && x>0 && x<1);
addParameter(p,'alpha_eq',0.95, @(x) isnumeric(x) && x>0 && x<1);
addParameter(p,'print',true, @islogical);
parse(p,varargin{:});
opts = p.Results;

% Clean chains
A = chainA(:); A = A(~isnan(A));
B = chainB(:); B = B(~isnan(B));
Na = length(A); Nb = length(B);
if Na==0 || Nb==0, error('Both chains must have at least one non-NaN sample.'); end

% Determine sampling strategy
if ~isempty(opts.sampling)
    sampling = char(opts.sampling);
else
    sampling = '';
end

% Determine nsamp default
if isempty(opts.nsamp)
    if Na == Nb
        nsamp = Na;        % use full paired set
    else
        nsamp = min(Na,Nb);% default to min length for deterministic pairing
    end
else
    nsamp = opts.nsamp;
end

% Build paired samples delta without randomness unless 'resample' chosen
switch sampling
    case {'paired'}
        if Na ~= Nb
            error('paired sampling requires equal-length chains.');
        end
        As = A;
        Bs = B;
        nsamp = Na;
        delta = As - Bs;
    case {'truncate'}
        m = min(Na,Nb);
        As = A(1:m);
        Bs = B(1:m);
        nsamp = m;
        delta = As - Bs;
    case {'resample'}
        % random resampling with replacement (user explicitly asked)
        if nsamp <= Na
            idxA = randsample(Na, nsamp, false);
        else
            idxA = randsample(Na, nsamp, true);
        end
        if nsamp <= Nb
            idxB = randsample(Nb, nsamp, false);
        else
            idxB = randsample(Nb, nsamp, true);
        end
        As = A(idxA);
        Bs = B(idxB);
        delta = As - Bs;
    otherwise % deterministic (default)
        if Na == Nb && nsamp == Na
            % elementwise paired
            As = A;
            Bs = B;
            delta = As - Bs;
        else
            % pick nsamp indices deterministically via linspace and round
            idxA = unique(round(linspace(1,Na,nsamp)));
            idxB = unique(round(linspace(1,Nb,nsamp)));
            % ensure equal length after uniqueness; if differ, adjust by truncation
            L = min(length(idxA), length(idxB));
            idxA = idxA(1:L);
            idxB = idxB(1:L);
            As = A(idxA);
            Bs = B(idxB);
            nsamp = L;
            delta = As - Bs;
        end
end

% Summary statistics
meanDiff = mean(delta);
lower = quantile(delta, (1-opts.cred)/2);
upper = quantile(delta, 1-(1-opts.cred)/2);
P_greater = mean(delta > 0);

% pooled SD
pooledSD = sqrt( ((Na-1)*var(A,0) + (Nb-1)*var(B,0)) / (Na+Nb-2) );

% ROPE
if strcmp(opts.mode,'std')
    if isempty(opts.rope)
        rope = 0.1 * pooledSD;
    else
        rope = opts.rope * pooledSD;
    end
else
    if isempty(opts.rope)
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

P_inROPE = mean(abs(delta) < rope);

% Decision
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
out.sampling = sampling;

% Print
if opts.print
    fprintf('Bayesian equality test (deterministic sampling)\n');
    fprintf(' sampling mode: %s\n', ternary(isempty(sampling),'auto',sampling));
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
    fprintf('\n\n')
end
end

% small helper: ternary
function out = ternary(cond, a, b)
    if cond, out = a; else out = b; end
end