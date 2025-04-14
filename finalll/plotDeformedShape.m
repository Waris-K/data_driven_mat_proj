% Function to plot deformed shape
function plotDeformedShape(nodes, elements, U, scaleFactor, binaryLayout)
    figure;
    hold on;
    
    % Plot original mesh
    for e = 1:size(elements, 1)
        nodeIndex = elements(e, :);
        patch(nodes(nodeIndex, 1), nodes(nodeIndex, 2), 'w', 'EdgeColor', 'b', 'FaceColor', 'none');
        
    end
    
    % Plot deformed mesh
    deformedNodes = nodes + scaleFactor * [U(1:2:end), U(2:2:end)];
    for e = 1:size(elements, 1)
        nodeIndex = elements(e, :);
        if binaryLayout(e)==1
            BlackorWhite = 'k';
        else
            BlackorWhite = 'w'; 
        end
        patch(deformedNodes(nodeIndex, 1), deformedNodes(nodeIndex, 2), 'w', 'EdgeColor', 'r', 'FaceColor', BlackorWhite);
      
    end
    
    title('Deformed Shape  vs Original Shape ');
    xlabel('X');
    ylabel('Y');
    axis equal;
    hold off;
end