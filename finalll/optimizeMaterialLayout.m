
function [x_opt, fval, strain_energy_fem, u_avg]=optimizeMaterialLayout()
% Uses a genetic algorithm to find an optimal layout

load('trained_surrogate.mat', 'net');

% Objective function
objective = @(x) predict(net, x);
%%
opts = optimoptions('ga','MaxStallGenerations',50,'FunctionTolerance',1e-10,...
    'PopulationSize',1000,'MaxGenerations',500,'PlotFcn',@gaplotbestfun);


% % Constraint: exactly 32 stiff elements
% nonlcon = @(x) deal([], sum(x) - 8);

% Binary variables
lb = zeros(16, 1);
ub = ones(16, 1);
intcon = 1:16;

rng(0,'twister') % for reproducibility
[x_opt,fval,exitflag] = ga(objective,16,[],[],[],[],lb,ub,@eqCon,intcon,opts);
% [x_opt, fval] = ga(objective, 16, [], [], [], [], lb, ub, nonlcon, intcon, options);

disp(['Best compliance found: ', num2str(fval)]);
[U, strain_energy_fem, u_avg] = runFEM1(x_opt);

% plotDeformedShape(U);
strain_energy_fem
u_avg
fval
x_opt'


% x_opt = [0     0     1     0     0     1     0     1     0     1     0     1     1     1     0     1];
grid_matrix = reshape(x_opt, 4, 4);
% grid_matrix = flipud(grid_matrix);
% grid_matrix = grid_matrix'; 
% grid_matrix = flipud(grid_matrix)
figure;
grid on
imagesc(grid_matrix);
colormap(gray)
axis equal tight;
end
% x_opt = [0     0     1     0     0     1     0     1     0     1     0     1     1     1     0     1];