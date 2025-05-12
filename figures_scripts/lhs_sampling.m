function samples = lhs_sampling(lower_lims, upper_lims, num_draws)

  % Check input dimensions
    if length(lower_lims) ~= length(upper_lims)
        error('lower_lims and upper_lims must be the same length');
    end
    
    N = length(lower_lims);  % number of parameters

    % Latin Hypercube Sampling in unit hypercube
    lhs_unit = lhsdesign(num_draws, N);  % num_draws x N matrix in [0, 1]

    % Scale samples to the given parameter bounds
    samples = bsxfun(@plus, lower_lims(:)', ...
              bsxfun(@times, lhs_unit, (upper_lims(:)' - lower_lims(:)')));
end
