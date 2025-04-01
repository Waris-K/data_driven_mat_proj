%% initialize
clear; clc; close all; 


VarName = 'material_distributionC.mat'; 
refinement_factor = 1; 

displacement = inplane_func(VarName, refinement_factor); 


save('displacement.mat',"displacement")

