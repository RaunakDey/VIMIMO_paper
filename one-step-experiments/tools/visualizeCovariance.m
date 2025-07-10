function visualizeCovariance(data,xlabels,ylabels)
    % Input: data - a 2-column matrix containing the sample data
    
    % Compute covariance matrix
    covMatrix = cov(data);
    
    % Extract the means of the two parameters
    mean1 = mean(data(:,1));
    mean2 = mean(data(:,2));
    
    % Create a grid for plotting
    [x, y] = meshgrid(linspace(min(data(:,1)), max(data(:,1)), 100), ...
                     linspace(min(data(:,2)), max(data(:,2)), 100));
    
    % Evaluate the bivariate normal PDF at each point on the grid
    z = mvnpdf([x(:) y(:)], [mean1, mean2], covMatrix);
    z = reshape(z, size(x));
    
    % Create the contour plot
    figure;
    contourf(x, y, z, 20); % Adjust the number of contours as needed
    colorbar;
     xlabel(string(xlabels));
    ylabel(string(ylabels));
    title('Covariance Contour Plot');
    set(gca,'FontSize',14)
    
    % Create the 3D plot
    figure;
    surf(x, y, z);
    colorbar;
    xlabel(xlabel);
    ylabel(ylabel);
    zlabel('Frequency');
    title('3D Covariance Plot');
end
