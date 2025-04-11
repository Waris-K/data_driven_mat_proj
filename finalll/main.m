
% Main script to generate dataset, train surrogate, and optimize layout
clear; clc; close all;
tic;
% Step 1: Generate dataset
numSamples = 50000;
disp('Generating FEM dataset...');
generateDataset(numSamples);
toc;
%%

% tic;
% Step 2: Train surrogate model
disp('Training surrogate model...');
trainSurrogate();
% toc;
%%


% tic;
% % 
% Step 3: Optimize layout
disp('Running optimization...');
[x_opt, fval, strain_energy_fem, u_avg]=optimizeMaterialLayout();
% toc;