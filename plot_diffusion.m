close all; clc; clear;
% Used to plot the 1D Velocity Field
[u_numerical, ~] = burgers_solve(150, 150, 0.5, 2.0*pi, 0.1);
[u_analytical, x] = analytical_solution(150, 150, 0.5, 2.0*pi, 0.1);

% Get each fourth term in the numeric solution so the plot is easy to see.
u_numeric_plot = u_numerical(1:4:end,:);
every_other_x = x(1:4:end);

% Initialize the legend data
legend_sub = [];
legend_labels = {};

% Iterate through each matrix and plot the numerical solution over
% the analytic solution at that iteration.
figure
for i = 1:5
    % Since n is large, extract the 40th term from the solution matricies.
    n = 30 * i;
    hold on
    
    % Store information for the legend
    legend_sub = [legend_sub plot(x, u_analytical(:, n))];
    legend_sub = [legend_sub scatter(every_other_x, u_numeric_plot(:, n), '*')];
    
    legend_labels{end + 1} = sprintf('Plot of the Analytic Solution After %i iterations', n);
    legend_labels{end + 1} = sprintf('Scatterplot of the Numeric Solution After %i iterations', n);

end
hold off

% Show the plot with the legend
set(gcf,'position',[600,500,1000,400])
legend(legend_sub, legend_labels);
legend('Location', 'northeastoutside');
xlabel('x (radians)')
ylabel('u (rad/s)');
title('Comparison of the Numerical and Analytical Solutions of the Burger PDE')

norms = zeros(1, length(u_analytical));

% Animate the solutions over time on the same plot
figure
hold on

% Store the numerical and analytical solutions separately
u_n = animatedline('Color', 'r');
u_a = animatedline('Color', 'b');
for i = 1:length(u_numerical)-1
    un_y = u_numerical(:,i);
    ua_y = u_analytical(:,i);
    norms(i) = norm(un_y - ua_y);
    
    % Add the points to the same graph
    addpoints(u_n, x, un_y)
    addpoints(u_a, x, ua_y)
    
    % Add the legends and labels
    legend('Evolution of Numerical Solution Over Time','Evolution of Analytical Solution Over Time');
    xlabel('x radians');
    ylabel('u (rad/s)');
    drawnow
    
    % Clear the points for the next iteration(s)
    clearpoints(u_a)
    clearpoints(u_n)
end
legend('off');
% Plot the numerical and analytical solutions for the last iteration
plot(x, u_numerical(:,length(u_numerical)), 'r')
plot(x, u_analytical(:,length(u_analytical)), 'b')
legend('Evolution of Numerical Solution Over Time','Evolution of Analytical Solution Over Time');

hold off

figure
% Plot the norm at each iteration
plot(1:length(norms), norms)
title('Norms at each Iteration')
xlabel('Iteration Number')
ylabel('Absolute Error (rad/s)')
