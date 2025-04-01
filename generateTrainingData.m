%% initialize
tic
clear; clc; close all; 

numSamples = 1000; % generate 1000 samples

% Geometry of the plate
Lx = 1.0;        % Length of the plate in x-direction (m)
Ly = 1.0;        % Length of the plate in y-direction (m)

% Mesh generation
nx = 8;         % Number of elements along x-axis
ny = 8;         % Number of elements along y-axis
[nodes, elements] = rectangularQuadMesh(Lx, Ly, nx, ny);

% Number of nodes and elements
numNodes = size(nodes, 1);
numElements = size(elements, 1);


num_material1 = numElements / 2; % Half assigned to material 1
num_material2 = numElements / 2; % Half assigned to material 2


for sample = 1:numSamples

    % Create an array with equal number of 0s and 1s, then shuffle
    material_list = [zeros(1, num_material1), ones(1, num_material2)];
    material_list = material_list(randperm(numElements));
    material_distribution(:,:,sample) =  reshape(material_list, ny, nx);

end

%Save material_distribution variable for later use feeding in to the FEM
save('material_distribution.mat', 'material_distribution');

% ploting one random distribution
figure;
imagesc(material_distribution(:, :, randi(numSamples)));
colormap([1 0 0; 0 0 1]); % Red for material 1, Blue for material 2
colorbar;
title('Example Material Distribution');
toc
%Elapsed time is 1.294244 seconds.