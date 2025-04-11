function [nodes, elements] = rectangularQuadMesh(Lx, Ly, nx, ny)
    % Generate node coordinates
    [x, y] = meshgrid(linspace(0, Lx, nx+1), linspace(0, Ly, ny+1));
    nodes = [x(:), y(:)];
    
    % Initialize element connectivity matrix
    elements = zeros(nx*ny, 4);
    
    % Assign nodes to elements in counter-clockwise order
    for i = 1:ny
        for j = 1:nx
            % Counter-clockwise node ordering:
            % Bottom-left -> Bottom-right -> Top-right -> Top-left
            node1 = (i-1)*(nx+1) + j;         % Bottom-left
            node2 = (i-1)*(nx+1) + j + 1;     % Bottom-right
            node3 = i*(nx+1) + j + 1;         % Top-right
            node4 = i*(nx+1) + j;             % Top-left
            
            elements((i-1)*nx + j, :) = [node1, node2, node3, node4];
        end
    end
end