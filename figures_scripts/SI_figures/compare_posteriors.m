function result = compare_posteriors(traceA, traceB, cred_level)
    % compare_posteriors: Compare two posterior traces
    % Inputs:
    %   traceA, traceB : vectors of posterior samples
    %   cred_level     : e.g. 0.95 (default)
    % Output:
    %   result struct with meanDiff, CrI, P_greater

    if nargin < 3
        cred_level = 0.95;
    end
    
    % Ensure column vectors
    traceA = traceA(:);
    traceB = traceB(:);
    
    % Match lengths by resampling the smaller one
    N = max(length(traceA), length(traceB));
    if length(traceA) < N
        traceA = randsample(traceA, N, true);
    end
    if length(traceB) < N
        traceB = randsample(traceB, N, true);
    end
    
    % Compute differences
    delta = traceA - traceB;
    
    % Posterior summaries
    meanDiff = mean(delta);
    CrI = quantile(delta, [(1-cred_level)/2, 1-(1-cred_level)/2]);
    P_greater = mean(delta > 0);
    
    % Store results
    result.meanDiff = meanDiff;
    result.CrI = CrI;
    result.P_greater = P_greater;
end