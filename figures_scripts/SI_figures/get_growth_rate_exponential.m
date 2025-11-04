function [mean_r, std_r, r_vals] = get_growth_rate_exponential(time, data)
% GET_GROWTH_RATE_EXPONENTIAL
%   [mean_r, std_r, r_vals] = get_growth_rate_exponential(time, data)
%   Fits r for each replicate from B(t)=B0*exp(r t) by linear fit to ln(B).
%   - time: vector length Nt
%   - data: Nt x R or R x Nt (replicates)
%   - mean_r: mean of replicate r values (omit NaN)
%   - std_r : std  of replicate r values (omit NaN)
%   - r_vals: R x 1 vector of replicate r (NaN if insufficient data)
%
%   Uses only positive finite points per replicate.

    time = time(:);
    Nt = numel(time);
    if Nt < 2
        error('time must have at least 2 points');
    end

    % orient data so rows = time, cols = replicates
    [d1,d2] = size(data);
    if d1 == Nt
        mat = data;        % Nt x R
    elseif d2 == Nt
        mat = data.';      % Nt x R
    else
        error('One dimension of data must match length(time).');
    end

    [Nt_check, R] = size(mat);
    if Nt_check ~= Nt
        error('Internal orientation failure');
    end

    r_vals = nan(R,1);

    % prepare plot
    figure;
    hold on;
    xlabel('Time');
    ylabel('ln(B)');
    title('Per-replicate ln(B) and exponential fits');
    tplot = linspace(min(time), max(time), 200);

    colors = lines(R);

    for rep = 1:R
        y = mat(:,rep);
        good = isfinite(y) & (y > 0);
        if sum(good) < 2
            % insufficient data for this replicate
            fprintf('rep %d: insufficient positive data -> NaN\n', rep);
            continue;
        end

        t_valid = time(good);
        lnB = log(y(good)); % natural log

        % linear fit lnB = a + r*t
        p = polyfit(t_valid, lnB, 1);
        r_vals(rep) = p(1);
        a = p(2);

        % plot data and fit
        plot(t_valid, lnB, 'o', 'Color', colors(rep,:), 'MarkerFaceColor', colors(rep,:), 'MarkerSize',6);
        plot(tplot, polyval(p, tplot), '-', 'Color', colors(rep,:), 'LineWidth', 1.25);

        % legend entry
        leg{rep} = sprintf('rep %d: r=%.4f', rep, r_vals(rep)); %#ok<AGROW>
    end

    % stats
    mean_r = mean(r_vals, 'omitnan');
    std_r  = std(r_vals, 'omitnan');

    % print summary
    if all(isnan(r_vals))
        fprintf('No valid replicates. mean_r and std_r = NaN\n');
    else
        fprintf('Per-replicate r (1/time):\n');
        for rep = 1:R
            if isnan(r_vals(rep))
                fprintf('  rep %d: NaN\n', rep);
            else
                fprintf('  rep %d: %.6f\n', rep, r_vals(rep));
            end
        end
        fprintf('Mean r = %.6f ± %.6f (mean ± std)\n', mean_r, std_r);
    end

    % finalize plot
    if exist('leg','var')
        legend(leg, 'Location', 'best');
    end
    % annotate mean±std
    ylimv = ylim;
    xpos = min(time) + 0.02*(max(time)-min(time));
    ypos = ylimv(2) - 0.06*(ylimv(2)-ylimv(1));
    text(xpos, ypos, sprintf('mean r = %.4f \\pm %.4f', mean_r, std_r), 'BackgroundColor','w','EdgeColor','k');

    hold off;
end