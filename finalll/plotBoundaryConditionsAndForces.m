% Function to plot boundary conditions and applied forces
function plotBoundaryConditionsAndForces(nodes, fixedNodes, forceNodes)
    figure;
    hold on;
    
    % Plot nodes
    plot(nodes(:, 1), nodes(:, 2), 'ko', 'MarkerSize', 5);
    
    % Plot fixed nodes (boundary conditions)
    plot(nodes(fixedNodes, 1), nodes(fixedNodes, 2), 'bs', 'MarkerSize', 10, 'LineWidth', 2);
    
    % Plot force nodes
    plot(nodes(forceNodes, 1), nodes(forceNodes, 2), 'r^', 'MarkerSize', 10, 'LineWidth', 2);
    
    % Add legend
    legend('Nodes', 'Fixed Nodes', 'Force Nodes');
    
    title('Boundary Conditions and Applied Forces');
    xlabel('X');
    ylabel('Y');
    axis equal;
    hold off;
end