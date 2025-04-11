% FEM code for in-plane loaded plate using quadrilateral elements (2D plane stress)
clear all; clc; close all;

% Material properties
E = 210e9;       % Young's modulus (Pa)
nu = 0.3;        % Poisson's ratio
t = 0.01;        % Thickness of the plate (m)

% Geometry of the plate
Lx = 1.0;        % Length of the plate in x-direction (m)
Ly = 1.0;        % Length of the plate in y-direction (m)

% Mesh generation
nx = 10;         % Number of elements along x-axis
ny = 10;         % Number of elements along y-axis
[nodes, elements] = rectangularQuadMesh(Lx, Ly, nx, ny);

% Number of nodes and elements
numNodes = size(nodes, 1);
numElements = size(elements, 1);

% Initialize global stiffness matrix and force vector
K = zeros(2*numNodes, 2*numNodes);
F = zeros(2*numNodes, 1);

% Gaussian quadrature points and weights
[gaussPoints, gaussWeights] = getGaussQuadrature(2);

% Assembly of global stiffness matrix
for e = 1:numElements
    % Element nodes
    nodeIndex = elements(e, :);
    elemNodes = nodes(nodeIndex, :);
    
    % Element stiffness matrix
    Ke = zeros(8, 8); % 4 nodes, 2 DOFs per node
    for gp = 1:length(gaussWeights)
        xi = gaussPoints(gp, 1);
        eta = gaussPoints(gp, 2);
        [N, dN_dxi, dN_deta] = shapeFunctions(xi, eta);
        J = Jacobian(elemNodes, dN_dxi, dN_deta);
        B = strainDisplacementMatrix(dN_dxi, dN_deta, J);
        D = materialMatrix(E, nu);
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

% Apply boundary conditions (fixed at left edge)
fixedNodes = find(nodes(:, 1) == 0); % Nodes on the left edge
fixedDofs = [2*fixedNodes-1; 2*fixedNodes]; % Corresponding DOFs
freeDofs = setdiff(1:2*numNodes, fixedDofs);

% Apply external force (uniform load on right edge)
forceNodes = find(nodes(:, 1) == Lx); % Nodes on the right edge
F(2*forceNodes) = 100000; % Apply force in (2*forceNodes-1) x-direction  (2*forceNodes) y-direction
F(2*forceNodes-1) = 100000;
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





