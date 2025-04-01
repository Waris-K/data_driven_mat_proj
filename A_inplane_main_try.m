% FEM code for in-plane loaded plate using quadrilateral elements (2D plane stress)
clear all; clc; close all;

%% load mat distribution
load('material_distribution.mat'); 

%% refine 
% Refinement Factor
refinementFactor = 1;  % Adjust to refine mesh (e.g., 2 = 16x16 elements)

% Compute New Mesh Size
size_mat = size(material_distribution); 
numElementsY = size_mat(1); 
numElementsX = numElementsY; 
newNumElementsX = numElementsX * refinementFactor;
newNumElementsY = numElementsY * refinementFactor;



% Generate Refined Material Distribution
new_material_distribution = zeros(newNumElementsY, newNumElementsX);

for row = 1:numElementsY
    for col = 1:numElementsX
        matType = material_distribution(row, col); % Get material type

        % Map to refined elements
        newRowStart = (row - 1) * refinementFactor + 1;
        newColStart = (col - 1) * refinementFactor + 1;

        % Assign the same material to all refined sub-elements
        new_material_distribution(newRowStart:newRowStart+refinementFactor-1, ...
                                  newColStart:newColStart+refinementFactor-1) = matType;
    end
end

%% Plot Refined Material Distribution
figure;
imagesc(new_material_distribution);
colormap(gray); % 0=dark, 1=light
title('Refined Material Distribution');
axis equal; axis off;


%%
% Material properties
E1 = 1000e6;       % Young's modulus (Pa)
E2 = 100e6;       % Young's modulus (Pa) % change back to 100

    % Note: this is an input to materialMatrix.m function need to adjust so
    % two moduli can be input for the two phases

nu = 0.3;        % Poisson's ratio
t = 0.01;        % Thickness of the plate (m)

% Geometry of the plate
Lx = 1.0;        % Length of the plate in x-direction (m)
Ly = 1.0;        % Length of the plate in y-direction (m)

% Mesh generation
nx = newNumElementsX;    % Number of elements along x-axis
ny = nx;         % Number of elements along y-axis
[nodes, elements] = rectangularQuadMesh(Lx, Ly, nx, ny);

% Number of nodes and elements
numNodes = size(nodes, 1);
numElements = size(elements, 1);

% Initialize global stiffness matrix and force vector
K = zeros(2*numNodes, 2*numNodes);
F = zeros(2*numNodes, 1);

% Gaussian quadrature points and weights
[gaussPoints, gaussWeights] = getGaussQuadrature(2);

%% Old Global Stiffness matrix generation and assembly
% % Assembly of global stiffness matrix
% for e = 1:numElements
%     % Element nodes
%     nodeIndex = elements(e, :);
%     elemNodes = nodes(nodeIndex, :);
%     
%     % Element stiffness matrix
%     Ke = zeros(8, 8); % 4 nodes, 2 DOFs per node
%     for gp = 1:length(gaussWeights)
%         xi = gaussPoints(gp, 1);
%         eta = gaussPoints(gp, 2);
%         [N, dN_dxi, dN_deta] = shapeFunctions(xi, eta);
%         J = Jacobian(elemNodes, dN_dxi, dN_deta);
%         B = strainDisplacementMatrix(dN_dxi, dN_deta, J);
%         D = materialMatrix(E, nu);
%         Ke = Ke + B' * D * B * det(J) * gaussWeights(gp) * t;
%     end
%     
%     % Assembly into global stiffness matrix
%     for i = 1:4
%         for j = 1:4
%             K(2*nodeIndex(i)-1:2*nodeIndex(i), 2*nodeIndex(j)-1:2*nodeIndex(j)) = ...
%                 K(2*nodeIndex(i)-1:2*nodeIndex(i), 2*nodeIndex(j)-1:2*nodeIndex(j)) + Ke(2*i-1:2*i, 2*j-1:2*j);
%         end
%     end
% end
%% Global Stiffness - changing E for different material 
% Assembly of global stiffness matrix


for e = 1:numElements
    % Element nodes
    nodeIndex = elements(e, :);
    elemNodes = nodes(nodeIndex, :);
    
    % Get material type for this element (using material_distribution matrix)
    % Assume material_distribution is an NxM matrix, where N=8, M=8 (size of the grid)
    [row, col] = ind2sub([newNumElementsY, newNumElementsX], e);
    material_type = new_material_distribution(row, col,1); % Extract material type from distribution matrix
    
    % Set the Young's Modulus based on the material type
    if material_type == 0
        E = E1; % Material 1
    else
        E = E2; % Material 2
    end
    
    % Element stiffness matrix
    Ke = zeros(8, 8); % 4 nodes, 2 DOFs per node
    for gp = 1:length(gaussWeights)
        xi = gaussPoints(gp, 1);
        eta = gaussPoints(gp, 2);
        [N, dN_dxi, dN_deta] = shapeFunctions(xi, eta);
        J = Jacobian(elemNodes, dN_dxi, dN_deta);
        B = strainDisplacementMatrix(dN_dxi, dN_deta, J);
        
        % Material matrix for current material (with E, nu)
        D = materialMatrix(E, nu);
        
        % Assemble element stiffness matrix
        Ke = Ke + B' * D * B * det(J) * gaussWeights(gp) * t;
    end
    
    % Assembly into global stiffness matrix
    for i = 1:4
        for j = 1:4
            K(2*nodeIndex(i)-1:2*nodeIndex(i), 2*nodeIndex(j)-1:2*nodeIndex(j)) = ...
                K(2*nodeIndex(i)-1:2*nodeIndex(i), 2*nodeIndex(j)-1:2*nodeIndex(j)) + Ke(2*i-1:2*i, 2*j-1:2*j);
        end
    end
end



%%


% Apply boundary conditions (fixed at left edge)
fixedNodes = find(nodes(:, 1) == 0); % Nodes on the left edge
fixedDofs = [2*fixedNodes-1; 2*fixedNodes]; % Corresponding DOFs
freeDofs = setdiff(1:2*numNodes, fixedDofs);

% % Apply external force (uniform load on right edge)
% forceNodes = find(nodes(:, 1) == Lx); % Nodes on the right edge
% % F(2*forceNodes) = 100000; % Apply force in (2*forceNodes-1) x-direction  (2*forceNodes) y-direction
% F(2*forceNodes-1) = -100;

% Apply external force (uniform load on right edge)
forceNodes = find(nodes(:, 1) == Lx); % Nodes on the right edge
% F(2*forceNodes) = 100000; % Apply force in (2*forceNodes-1) x-direction  (2*forceNodes) y-direction
F(2*forceNodes-1) = -1000/length(forceNodes); %distribute load based on number of nodes

% Solve for displacements
U = zeros(2*numNodes, 1);
U(freeDofs) = K(freeDofs, freeDofs) \ F(freeDofs);

% Calculate stresses
stresses = zeros(numElements, 3); % [sigma_xx, sigma_yy, tau_xy]
for e = 1:numElements
    nodeIndex = elements(e, :);
    elemNodes = nodes(nodeIndex, :);
    Ue = zeros(8, 1); % Initialize Ue as an 8x1 vector
    for i = 1:4
        Ue(2*i-1) = U(2*nodeIndex(i)-1); % x-displacement for node i
        Ue(2*i) = U(2*nodeIndex(i));     % y-displacement for node i
    end
    
    % Evaluate shape functions and derivatives at the centroid (xi = 0, eta = 0)
    [N, dN_dxi, dN_deta] = shapeFunctions(0, 0);
    
    % Compute Jacobian matrix
    J = Jacobian(elemNodes, dN_dxi, dN_deta);
    
    % Compute strain-displacement matrix B
    B = strainDisplacementMatrix(dN_dxi, dN_deta, J); % B should be 3x8
    
    % Compute material matrix D
    D = materialMatrix(E, nu); % D should be 3x3
    
    % Compute stresses: sigma = D * B * Ue
    stresses(e, :) = (D * B * Ue)'; 
end

% Display results
% disp('Nodal Displacements:');
% disp(U);
% disp('Element Stresses:');
% disp(stresses);
%%
% Plot boundary conditions and applied forces
plotBoundaryConditionsAndForces(nodes, fixedNodes, forceNodes);

% Plot deformed shape with scaling
scaleFactor = 100; % Scale factor for visualization
plotDeformedShape(nodes, elements, U, scaleFactor);

%% Reordering displacement and  plotting


% Number of nodes in each direction
numNodesX = nx + 1; 
numNodesY = ny + 1; 

% Assume 'U' is the displacement column vector (size: 2*numNodes x 1)
Ux = U(1:2:end); % Extract x-displacements
Uy = U(2:2:end); % Extract y-displacements

% Reshape into mesh-sized matrices
Ux_matrix = reshape(Ux, [numNodesY, numNodesX]); % ny rows, nx columns
Uy_matrix = reshape(Uy, [numNodesY, numNodesX]); % ny rows, nx columns

% Plot results (optional)
figure;
subplot(1,2,1);
imagesc(Ux_matrix);
colorbar;
title('X Displacement');

subplot(1,2,2);
imagesc(Uy_matrix);
colorbar;
title('Y Displacement');


