% Function to compute material matrix
function D = materialMatrix(E, nu)
    D = E / (1 - nu^2) * [1, nu, 0;
                          nu, 1, 0;
                          0, 0, (1 - nu)/2];
end