%% MATLAB 3D Plotting Script with Regression Analysis (Z vs Y)


% 2. Import the data
data = readmatrix('testxyz.csv');

% Extract X, Y, Z (Assuming cols 1, 2, 3)
x = data(:, 1);
y = data(:, 2);
z = data(:, 3);

% 3. Calculate Linear Regression: Z = f(Y)
% polyfit(independent, dependent, degree)
coeffs = polyfit(y, z, 1); 
slope = coeffs(1);
intercept = coeffs(2);

% Calculate R-squared (Correlation Coefficient)
R_matrix = corrcoef(y, z);
r_sq = R_matrix(1,2)^2;

% 4. Create the Visualization
figure('Color', 'w');
scatter3(x, y, z, 30, z, 'filled'); 
% 1. Calculate the regression line points
p = polyfit(y, z, 1);                  % Linear fit: z = p(1)*y + p(2)
y_line = linspace(min(y), max(y), 100); % Create range of Y values
z_line = polyval(p, y_line);           % Calculate corresponding Z values
%x_line = repmat(mean(x), size(y_line)); % Position line at the center of X
x_line = repmat(max(x), size(y_line));
% 2. Plot the line
hold on;                               % Keep the original data on the plot
plot3(x_line, y_line, z_line, 'r-', 'LineWidth', 3); 
hold off;
hold on; % Keep the scatter plot while adding regression info

% 5. Add Labels and Styling
grid on;
colormap(parula);
colorbar;
xlabel('Degree');
ylabel('Tau (max coeff bitsize)');
zlabel('Execution time');
title(['3D Plot: ', file], 'Interpreter', 'none');

% 6. Display Regression Results in a Text Box
statsStr = { ...
   ['Regression: Z = f(Y)'], ...
   ['Slope: ', num2str(slope, '%.4f')], ...
   ['R^2: ', num2str(r_sq, '%.4f')] ...
};

% Places text in the top-left of the figure (normalized units)
annotation('textbox', [0.15, 0.75, 0.2, 0.15], ...
   'String', statsStr, ...
   'BackgroundColor', 'white', ...
   'EdgeColor', 'black', ...
   'LineWidth', 1);

view(3);
axis tight;
hold off;
