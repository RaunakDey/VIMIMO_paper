function h = violinNchains(chains, varargin)
% VIOLINNCHAINS  Side-by-side violin plots for N chains
% h = violinNchains(chains, Name,Value...)
%
% Inputs
%  chains        : cell array {c1,c2,...} of vectors OR numeric matrix
%                  with columns as chains (each column is a chain)
% Name-Value options
%  'labels'      : cell array of labels (default: {'1','2',...})
%  'colors'      : Nx3 RGB matrix or single RGB row (default gray palette)
%  'width'       : max half-width of violins (default 0.3)
%  'bandwidth'   : bandwidth for ksdensity (default [])
%  'alpha'       : face alpha (default 0.7)
%  'showScatter' : true/false to plot samples (default false)
%  'support'     : ksdensity support (default [0 Inf]; use [] for unrestricted)
%
% Output
%  h : struct with fields .patch (handles), .median (handles), .scatter (handles)

% ---- parse inputs
if ismatrix(chains) && ~iscell(chains)
    % numeric matrix: treat columns as chains
    C = cell(1,size(chains,2));
    for k=1:size(chains,2), C{k} = chains(:,k); end
else
    C = chains;
end
N = numel(C);

p = inputParser;
addParameter(p,'labels', arrayfun(@(k) sprintf('%d',k), 1:N, 'UniformOutput', false));
addParameter(p,'colors', []);
addParameter(p,'width', 0.3, @(x) isnumeric(x) && isscalar(x));
addParameter(p,'bandwidth', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p,'alpha', 0.7, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
addParameter(p,'showScatter', false, @islogical);
addParameter(p,'support', [0 Inf]); % default nonnegative support
parse(p, varargin{:});
opts = p.Results;

% default colors (if not provided)
if isempty(opts.colors)
    base = repmat([0.5 0.5 0.5], N, 1);
    % add slight variation for subsequent chains
    for k=1:N
        base(k,:) = base(k,:) + 0.12 * (mod(k-1,3)-1);
    end
    colors = min(max(base,0),1);
else
    colors = opts.colors;
    if size(colors,1) == 1
        colors = repmat(colors, N, 1);
    elseif size(colors,1) < N
        % repeat rows cyclically
        colors = repmat(colors, ceil(N/size(colors,1)), 1);
        colors = colors(1:N,:);
    end
end

hold on
patchHandles = gobjects(N,1);
medianHandles = gobjects(N,1);
scatterHandles = gobjects(N,1);

for i = 1:N
    data = C{i}(:);
    data = data(~isnan(data));
    if isempty(data)
        warning('Chain %d is empty. Skipping.', i);
        continue
    end

    % density estimate
    if isempty(opts.bandwidth)
        if isempty(opts.support)
            [f,xi] = ksdensity(data);
        else
            [f,xi] = ksdensity(data,'Support',opts.support);
        end
    else
        if isempty(opts.support)
            [f,xi] = ksdensity(data,'Bandwidth',opts.bandwidth);
        else
            [f,xi] = ksdensity(data,'Bandwidth',opts.bandwidth,'Support',opts.support);
        end
    end

    % normalize to requested width
    if isempty(f) || max(f)==0
        f = zeros(size(xi));
    else
        f = f / max(f) * opts.width;
    end

    % polygon coordinates
    xpoly = [i - f, fliplr(i + f)];
    ypoly = [xi, fliplr(xi)];

    % draw patch
    patchHandles(i) = patch(xpoly, ypoly, colors(i,:), 'FaceAlpha', opts.alpha, 'EdgeColor', 'none');

    % median line
    med = median(data);
    medianHandles(i) = plot([i-opts.width, i+opts.width], [med med], 'k-', 'LineWidth', 1.5);

    % optional scatter of raw points (jittered)
    if opts.showScatter
        jitter = (rand(size(data)) - 0.5) * opts.width * 0.5;
        scatterHandles(i) = plot(i + jitter, data, '.', 'Color', colors(i,:)*0.6, 'MarkerSize', 6);
    else
        scatterHandles(i) = gobjects(1);
    end
end

xlim([0.5, N+0.5]);
set(gca, 'XTick', 1:N, 'XTickLabel', opts.labels);
box on
hold off

h.patch = patchHandles;
h.median = medianHandles;
h.scatter = scatterHandles;
end