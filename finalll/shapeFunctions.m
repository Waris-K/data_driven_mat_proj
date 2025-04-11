% Function to compute shape functions and derivatives
function [N, dN_dxi, dN_deta] = shapeFunctions(xi, eta)
    % Shape functions for 4-node quadrilateral element
    N = 0.25 * [
        (1 - xi) * (1 - eta);
        (1 + xi) * (1 - eta);
        (1 + xi) * (1 + eta);
        (1 - xi) * (1 + eta)
    ];
    
    % Derivatives of shape functions w.r.t. xi and eta
    dN_dxi = 0.25 * [
        -(1 - eta),  (1 - eta),  (1 + eta), -(1 + eta)
    ];
    dN_deta = 0.25 * [
        -(1 - xi), -(1 + xi),  (1 + xi),  (1 - xi)
    ];
end

