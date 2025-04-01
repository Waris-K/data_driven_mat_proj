clc; clear; close all;

% Grid size
Nx = 8; % Number of elements in X
Ny = 8; % Number of elements in Y
Lx = 1; % Total length in X
Ly = 1; % Total length in Y
dx = Lx / Nx; % Element width
dy = Ly / Ny; % Element height

% Define material properties
E1 = 200e9; % Young's modulus for material 1 (e.g., steel)
E2 = 100e9; % Young's modulus for material 2 (e.g., aluminum)

% Total number of elements
num_elements = Nx * Ny;
num_material1 = num_elements / 2; % Half assigned to material 1
num_material2 = num_elements / 2; % Half assigned to material 2

% Create an array with equal number of 0s and 1s, then shuffle
material_list = [zeros(1, num_material1), ones(1, num_material2)];
material_list = material_list(randperm(num_elements));

% Reshape into grid form
material_assignment = reshape(material_list, Ny, Nx);



% Generate mesh grid
x = linspace(0, Lx, Nx+1);
y = linspace(0, Ly, Ny+1);
[X, Y] = meshgrid(x, y);

%iterate material property
E = E1 * (material_assignment == 0) + E2 * (material_assignment == 1);

% Plot the undeformed shape
figure;
hold on;
for i = 1:Nx
    for j = 1:Ny
        % Get element corners
        x_elem = [x(i), x(i+1), x(i+1), x(i)];
        y_elem = [y(j), y(j), y(j+1), y(j+1)];
        
        % Assign color based on material
        if material_assignment(j, i) == 0
            fill(x_elem, y_elem, 'r'); % Material 1 (Red)
        else
            fill(x_elem, y_elem, 'b'); % Material 2 (Blue)
        end
    end
end

xlabel('X'); ylabel('Y');
title('Undeformed Mesh with Random Material Distribution');
axis equal;
grid on;
hold off;
