function y = gaussian_truncated(x, mu, sigma, lower, upper)
    % Check inputs
    if sigma <= 0
        error('Standard deviation sigma must be positive.');
    end
    if lower >= upper
        error('Lower bound must be less than upper bound.');
    end
    
    % Unnormalized PDF
    y_raw = 1/(sqrt(2*pi)*sigma) .* exp(-0.5*((x - mu)./sigma).^2);
    
    % Normalization constant using the CDF
    Z = normcdf(upper, mu, sigma) - normcdf(lower, mu, sigma);
    
    % Truncate: Set PDF to 0 outside [lower, upper]
    y = y_raw ./ Z;
    y(x < lower | x > upper) = 0;
end