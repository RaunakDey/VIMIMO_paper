function violin2chains(chainA, chainB, varargin)
% VIOLIN2CHAINS  Plot two chains as side-by-side violins
% violin2chains(chainA, chainB, 'labels', {'Model A','Model B'}, 'colors', [0 0.45 0.74; 0.85 0.33 0.1])
%
% Inputs:
%   chainA, chainB : vectors of samples
% Optional (name-value):
%   'labels'  : cell array of two strings
%   'colors'  : 2x3 matrix of RGB values
%   'width'   : maximum violin half-width (default 0.3)
%   'bandwidth': kernel bandwidth for ksdensity (default auto)

% Parse inputs
p = inputParser;
addParameter(p,'labels',{'Chain A','Chain B'});
%addParameter(p,'colors',[0.2 0.5 0.8; 0.8 0.3 0.3]);
addParameter(p,'colors',[0.5 0.5 0.5; 171/255 193/255 157/255]);

addParameter(p,'width',0.3);
addParameter(p,'bandwidth',[]);
parse(p,varargin{:});
opts = p.Results;

chains = {chainA(:), chainB(:)};
hold on
for i = 1:2
    data = chains{i};
    data = data(~isnan(data));
    %if isempty(opts.bandwidth)
    %    [f,xi] = ksdensity(data);
    %else
    %    [f,xi] = ksdensity(data,'bandwidth',opts.bandwidth);
    %end
    
    % compute density restricted to [0, Inf)
    if isempty(opts.bandwidth)
        [f,xi] = ksdensity(data,'Support',[0 Inf]);
    else
        [f,xi] = ksdensity(data,'Support',[0 Inf],'Bandwidth',opts.bandwidth);
    end

    f = f / max(f) * opts.width; % normalize to width
    
    % Coordinates for violin
    x = [i - f, fliplr(i + f)];
    y = [xi, fliplr(xi)];
    
    % Patch (violin)
    patch(x, y, opts.colors(i,:), 'FaceAlpha',0.7, 'EdgeColor','none');
    
    % Median line
    med = median(data);
    plot([i-opts.width, i+opts.width], [med med], 'k-', 'LineWidth',1.5);
    
    % Scatter raw samples (jittered horizontally)
    jitter = (rand(size(data))-0.5) * opts.width*0.5;
    %plot(i + jitter, data, '.', 'Color', opts.colors(i,:)*0.6, 'MarkerSize',4);
end

xlim([0.5 2.5])
set(gca,'XTick',[])%,'XTickLabel',opts.labels)
%ylabel('Value')
%title('Violin plot comparison')
box on
hold off
end