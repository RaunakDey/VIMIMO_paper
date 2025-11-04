function result = compare_posteriors_equiv(traceA, traceB, cred_level, rope)
    % Compare two posteriors for equality vs difference
    %
    % Inputs:
    %   traceA, traceB : posterior sample vectors
    %   cred_level     : e.g. 0.95
    %   rope           : tolerance for practical equivalence
    %
    % Output: struct with meanDiff, CrI, P_inROPE

    if nargin < 3
        cred_level = 0.95;
    end
    if nargin < 4
        rope = 0; % default: exact equality
    end
    
    traceA = traceA(:); traceB = traceB(:);
    N = max(length(traceA), length(traceB));
    if length(traceA) < N, traceA = randsample(traceA, N, true); end
    if length(traceB) < N, traceB = randsample(traceB, N, true); end
    
    delta = traceA - traceB;
    
    meanDiff = mean(delta);
    CrI = quantile(delta, [(1-cred_level)/2, 1-(1-cred_level)/2]);
    P_inROPE = mean(abs(delta) < rope);
    
    result.meanDiff = meanDiff;
    result.CrI = CrI;
    result.P_inROPE = P_inROPE;
end