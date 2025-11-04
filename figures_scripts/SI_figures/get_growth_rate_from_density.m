function [slope_log10, slope_ln] = get_growth_rate_from_density(time, data)
% GET_GROWTH_RATE_FROM_DENSITY  Fit exponential growth to density data.
%   [slope_log10, slope_ln] = get_growth_rate_from_density(time, data)
%   - time: vector length Nt (can be row or column).
%   - data: Nt x R or R x Nt matrix (Nt = number of timepoints,
%           R = number of replicates). Function autodetects orientation.
%   - slope_log10: fitted slope in units log10(density)/time.
%   - slope_ln: slope in natural-log units (1/time) = slope_log10 * ln(10).
%
%   Behavior:
%   - At each timepoint take mean across replicates using only finite >0 values.
%   - Fit linear model to log10(meanDensity) vs time using polyfit.
%   - Returns NaN if fewer than 2 valid timepoints.

    % normalize time
    time = time(:);
    Nt = numel(time);
    if Nt < 2
        error('time must have at least 2 points');
    end

    % orient data so rows = timepoints, cols = replicates
    [d1,d2] = size(data);
    if d1 == Nt
        mat = data;            % Nt x R
    elseif d2 == Nt
        mat = data.';         % Nt x R
    else
        error('One dimension of data must match length(time).');
    end

    % compute mean across replicates using only positive finite values
    meanDensity = nan(Nt,1);
    for i = 1:Nt
        vals = mat(i,:);
        good = isfinite(vals) & (vals > 0);
        if any(good)
            meanDensity(i) = mean(vals(good));
        else
            meanDensity(i) = NaN;
        end
    end

    valid = isfinite(meanDensity) & (meanDensity > 0);
    if sum(valid) < 2
        slope_log10 = NaN;
        slope_ln = NaN;
        fprintf('Insufficient valid averaged points. Returning NaN.\n');
        return;
    end

    t_valid = time(valid);
    ylog = log10(meanDensity(valid));

    % linear fit in log10 space
    p = polyfit(t_valid, ylog, 1);
    slope_log10 = p(1);
    slope_ln = slope_log10 * log(10); % convert to natural-log units (1/time)

    % print
    fprintf('Growth rate = %.6f (log10 units/time) = %.6f (1/time, ln units)\n', ...
            slope_log10, slope_ln);

    % plot
    figure;
    hold on;
    plot(t_valid, ylog, 'ko', 'MarkerFaceColor','k');
    tplot = linspace(min(t_valid), max(t_valid), 200);
    plot(tplot, polyval(p, tplot), 'LineWidth', 1.25);
    xlabel('Time');
    ylabel('log_{10}(mean density)');
    title('Averaged log_{10}(density) and linear fit');
    txt = sprintf('slope_{log10}=%.4f, slope_{ln}=%.4f', slope_log10, slope_ln);
    ylimv = ylim;
    text(min(t_valid)+0.02*(max(t_valid)-min(t_valid)), ...
         ylimv(2)-0.06*(ylimv(2)-ylimv(1)), txt, ...
         'BackgroundColor','w','EdgeColor','k');
    hold off;
end