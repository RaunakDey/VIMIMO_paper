% Define the binary matrix
binary_matrix = [ 0 1 0 0 0;
                 1 1 1 0 0;
                 0 0 1 0 0;
                 0 0 0 1 1;
                 0 0 0 1 1];

% Define the color for value '1'
color_blue = [0.5 0.5 0.5];

% Create the figure
figure;
hold on;

count = 0;
% Loop through each element in the matrix
for j = 1:5
    for i = 1:5
        % Determine the color based on the value
        if binary_matrix(i, j) == 1
            face_color = color_blue;
            count = count + 1; 
        else
            face_color = [1 1 1]; % White for '0'
        end
        

     
        % Draw the rounded rectangle
        rectangle('Position', [j-0.9, 5-i+0.1, 0.8, 0.8], ...
                  'Curvature', [0.2, 0.2], ...
                  'EdgeColor', 'k', ... % Black boundary
                  'FaceColor', face_color, ...
                  'LineWidth', 1.5);
        
        if binary_matrix(i, j) == 1
        % Add the text inside the box
        text(j-0.5, 5-i+0.5, num2str(count), ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle','FontSize',36,'FontName','times',Color=[1 1 1]);
        end

    end
end

% Set x-axis and y-axis labels


set(gca, 'XTick', (1:5) -0.5, 'XTickLabel', {'\phi18:2', '\phi18:3', '\phi38:1', 'PSA-HP1', 'PSA-HS6'});
xtickangle(90);

set(gca, 'YTick', (1:5) -0.5, 'YTickLabel', {'PSA 13-15','PSA H100','CBA 38', 'CBA 18','CBA 4'});




 set(gca,'Fontsize',28);
set(gca,'fontname','times');

% Set the axis limits and reverse the y-axis direction
axis([0 5 0 5]);
axis("square")
%set(gca, 'YDir', 'reverse');
set(gca,'xaxisLocation','top')


% Add grid lines

% Turn off the hold
hold off;
box off;

%%

ax = gca;
% Set color of X and Y axes to background color (white)
set(ax,'XColor',[1 1 1])
set(ax,'YColor',[1 1 1])
% Set color of axis labels back to standard color
set(ax.XLabel,'Color',[0.15 0.15 0.15])
set(ax.YLabel,'Color',[0.15 0.15 0.15])
% Set color of tick labels back to standard color
set(ax.XAxis,'TickLabelColor',[0.15 0.15 0.15])
set(ax.YAxis,'TickLabelColor',[0.15 0.15 0.15])

