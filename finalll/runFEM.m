
function [U, strain_energy] = runFEM(materialLayout)
% FEM simulation for a 1x1 plate with 64 square quads and 2 materials

E_stiff = 1e9;
E_soft = 0.1e9;
nu = 0.3;
t = 0.01;
Lx = 1.0; Ly = 1.0;
nx = 2; ny = 2;

[nodes, elements] = rectangularQuadMesh(Lx, Ly, nx, ny);
numNodes = size(nodes, 1);
numElements = size(elements, 1);
K = zeros(2*numNodes);
F = zeros(2*numNodes,1);
[gp, wts] = getGaussQuadrature(2);

for e = 1:numElements
    nodeIdx = elements(e,:);
    elemNodes = nodes(nodeIdx, :);
    Ke = zeros(8);
    for gpId = 1:length(wts)
        xi = gp(gpId,1); eta = gp(gpId,2);
        [~, dNxi, dNeta] = shapeFunctions(xi, eta);
        J = Jacobian(elemNodes, dNxi, dNeta);
        B = strainDisplacementMatrix(dNxi, dNeta, J);
        mat = materialMatrix((materialLayout(e)==1)*E_stiff + (materialLayout(e)==0)*E_soft, nu);
        Ke = Ke + B'*mat*B*det(J)*wts(gpId)*t;
    end
    for i = 1:4
        for j = 1:4
            K(2*nodeIdx(i)-1:2*nodeIdx(i), 2*nodeIdx(j)-1:2*nodeIdx(j)) = ...
                K(2*nodeIdx(i)-1:2*nodeIdx(i), 2*nodeIdx(j)-1:2*nodeIdx(j)) + Ke(2*i-1:2*i, 2*j-1:2*j)
        end
    end
end

fixedNodes = find(nodes(:,1)==0);
fixedDofs = [2*fixedNodes-1; 2*fixedNodes];
freeDofs = setdiff(1:2*numNodes, fixedDofs);

forceNodes = find(nodes(:,1)==Lx);
F(2*forceNodes-1) = 1e5;

U = zeros(2*numNodes,1);
U(freeDofs) = K(freeDofs, freeDofs) \ F(freeDofs);



strain_energy = 0.5 *U'*K*U


end
