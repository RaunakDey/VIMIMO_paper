function growth_rate = exponential_fit_averaged(time, OD)
% EXPONENTIAL_FIT_AVERAGED
%   growth_rate = exponential_fit_averaged(time, OD)
%   - time: vector length N
%   - OD: R x N or N x R matrix (R replicates)
%   - Averages positive finite OD across replicates at each timepoint,
%     fits log10(meanOD) = a + b*t, returns b (log10 units per time).
%
%   If insufficient valid points, returns NaN.

    % normalize time
    time = time(:);
    if numel(time) < 2
        error('time must have at least 2 points');
    end

    % normalize OD to R x N
    [r1, c1] = size(OD);
    if r1 == numel(time) && c1 ~= numel(time)
        OD = OD';
    end
    [R, N] = size(OD);
    if N ~= numel(time)
        error('length(time) must match number of columns in OD');
    end

    % compute mean across replicates at each time using only positive finite values
    meanOD = nan(1, N);
    for j = 1:N
        vals = OD(:, j);
        good = isfinite(vals) & (vals > 0);
        if any(good)
            meanOD(j) = mean(vals(good));
        else
            meanOD(j) = NaN;
        end
    end

    % select timepoints with valid meanOD
    valid = isfinite(meanOD) & (meanOD > 0);
    if sum(valid) < 2
        growth_rate = NaN;
        fprintf('Insufficient valid averaged OD points. Output = NaN\n');
        return;
    end

    t_valid = time(valid);
    logy = log10(meanOD(valid));

    % linear fit in log10 space
    p = polyfit(t_valid, logy, 1);
    b = p(1); a = p(2);
    growth_rate = b; % units: log10(OD) per time

    % print result
    fprintf('Growth rate (log10 units/time) = %.6f\n', growth_rate);

    % plot
    figure;
    hold on;
    plot(t_valid, logy, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    tplot = linspace(min(t_valid), max(t_valid), 200)';
    plot(tplot, polyval(p, tplot), 'LineWidth', 1.25);
    xlabel('Time');
    ylabel('log_{10}(mean OD)');
    title('Averaged log_{10}(OD) and linear fit');
    txt = sprintf('growth rate = %.4f (log10 units/time)', growth_rate);
    ylim_vals = ylim;
    text_x = min(t_valid) + 0.02*(max(t_valid)-min(t_valid));
    text_y = ylim_vals(2) - 0.06*(ylim_vals(2)-ylim_vals(1));
    text(text_x, text_y, txt, 'BackgroundColor', 'w', 'EdgeColor', 'k');
    hold off;
end