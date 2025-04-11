function J = Jacobian(elemNodes, dN_dxi, dN_deta)
    % Compute Jacobian matrix (2x2)
    J = zeros(2, 2);
    J(1,1) = dN_dxi * elemNodes(:,1);   % ∂x/∂ξ
    J(1,2) = dN_dxi * elemNodes(:,2);   % ∂y/∂ξ
    J(2,1) = dN_deta * elemNodes(:,1);  % ∂x/∂η
    J(2,2) = dN_deta * elemNodes(:,2);  % ∂y/∂η
end