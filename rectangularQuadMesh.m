% Function to generate rectangular mesh with quadrilateral elements
function [nodes, elements] = rectangularQuadMesh(Lx, Ly, nx, ny)
    [x, y] = meshgrid(linspace(0, Lx, nx+1), linspace(0, Ly, ny+1));
    nodes = [x(:), y(:)];
    elements = zeros(nx*ny, 4);
    for i = 1:ny
        for j = 1:nx
            elements((i-1)*nx + j, :) = [
                (i-1)*(nx+1) + j, ...
                (i-1)*(nx+1) + j + 1, ...
                i*(nx+1) + j + 1, ...
                i*(nx+1) + j
            ];
        end
    end
end
