% Function to get Gaussian quadrature points and weights
function [points, weights] = getGaussQuadrature(order)
    if order == 1
        points = [0, 0];
        weights = 4;
    elseif order == 2
        points = [-1/sqrt(3), -1/sqrt(3);
                   1/sqrt(3), -1/sqrt(3);
                  -1/sqrt(3),  1/sqrt(3);
                   1/sqrt(3),  1/sqrt(3)];
        weights = [1; 1; 1; 1];
    end
end