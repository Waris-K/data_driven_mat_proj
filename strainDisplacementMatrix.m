% Function to compute strain-displacement matrix
function B = strainDisplacementMatrix(dN_dxi, dN_deta, J)
    % Compute derivatives of shape functions w.r.t. global coordinates
    dN_dx = inv(J) * [dN_dxi; dN_deta]; % 2x4 matrix
    
    % Initialize B matrix (3x8)
    B = zeros(3, 8);
    
    % Fill B matrix
    for i = 1:4
        B(1, 2*i-1) = dN_dx(1, i); % epsilon_xx
        B(2, 2*i)   = dN_dx(2, i); % epsilon_yy
        B(3, 2*i-1) = dN_dx(2, i); % gamma_xy
        B(3, 2*i)   = dN_dx(1, i); % gamma_xy
    end
end