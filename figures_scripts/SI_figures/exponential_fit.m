function [mean_rate, std_rate, growth_rates] = exponential_fit(time, OD)
% EXPONENTIAL_FIT Fit exponential growth to OD600 replicates (log10 scale).
%   [mean_rate, std_rate, growth_rates] = exponential_fit(time, OD)
%   - time: vector length N (e.g. 28)
%   - OD:  R x N or N x R matrix (R replicates)
%   - returns mean and std of fitted slopes (growth rates, units = log10(OD)/time)
%   - also returns growth_rates (R x 1), prints mean±std and per-replicate rates
%   - plots log10(OD) and fitted lines

    % normalize inputs
    time = time(:);
    if numel(time) < 2
        error('time must have at least 2 points');
    end

    [r1, c1] = size(OD);
    if r1 == numel(time) && c1 ~= numel(time)
        % OD is N x R, transpose to R x N
        OD = OD';
    end

    [R, N] = size(OD);
    if N ~= numel(time)
        error('length(time) must match number of columns in OD');
    end

    growth_rates = nan(R,1);
    intercepts = nan(R,1);

    % Prepare plot
    figure;
    hold on;
    xlabel('Time');
    ylabel('log_{10}(OD)');
    title('Log_{10}-transformed OD and linear fits');
    legends = cell(0);

    tplot = linspace(min(time), max(time), 200)';

    for r = 1:R
        y = OD(r, :)';
        valid = isfinite(y) & (y > 0);    % require positive OD for log10
        if sum(valid) < 2
            continue;
        end

        t_valid = time(valid);
        logy = log10(y(valid));

        % Linear fit: log10(OD) = a + b*t
        p = polyfit(t_valid, logy, 1);
        b = p(1); a = p(2);
        growth_rates(r) = b;
        intercepts(r) = a;

        % Plot points and fitted line
        plot(t_valid, logy, 'o', 'MarkerSize', 6);
        plot(tplot, polyval(p, tplot), 'LineWidth', 1.25);
        legends{end+1} = sprintf('rep %d: %.4f', r, b); %#ok<AGROW>
    end

    % Compute stats
    mean_rate = mean(growth_rates, 'omitnan');
    std_rate  = std(growth_rates, 'omitnan');

    % Print summary
    if all(isnan(growth_rates))
        fprintf('No valid replicates with positive OD found. Outputs are NaN.\n');
    else
        fprintf('Per-replicate growth rates (log10 units / time):\n');
        for r = 1:R
            if isnan(growth_rates(r))
                fprintf('  rep %d: insufficient positive data\n', r);
            else
                fprintf('  rep %d: %.6f\n', r, growth_rates(r));
            end
        end
        fprintf('Mean growth rate = %.6f ± %.6f (mean ± std)\n', mean_rate, std_rate);
    end

    % Annotate plot with mean±std
    text_x = min(time) + 0.02*(max(time)-min(time));
    text_y = max(ylim) - 0.06*(max(ylim)-min(ylim));
    txt = sprintf('mean = %.4f \\pm %.4f', mean_rate, std_rate);
    text(text_x, text_y, txt, 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k');

    if ~isempty(legends)
        legend(legends, 'Location', 'best');
    end
    hold off;
end