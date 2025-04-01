% Function to compute Jacobian matrix
function J = Jacobian(elemNodes, dN_dxi, dN_deta)
    J = [dN_dxi * elemNodes(:, 1), dN_dxi * elemNodes(:, 2);
         dN_deta * elemNodes(:, 1), dN_deta * elemNodes(:, 2)];
end