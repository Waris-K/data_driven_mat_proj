
function generateDataset(numSamples)
% Generates stiffness/compliance data for random layouts

data = [];
for i = 1:numSamples
    layout = randperm(16); % Random layout permutation
    binaryLayout = zeros(16,1);
    binaryLayout(layout(1:8)) = 1; % 1 for stiff, 0 for soft
    
    [~, strain_energy, u_avg] = runFEM1(binaryLayout);
    % data = [data; binaryLayout',U, strain_energy];
    data = [data; binaryLayout', strain_energy];
    % data = [data; binaryLayout', u_avg];
end

save('fem_dataset.mat', 'data');
end

% X = zeros(numSamples, 64);
% Y = zeros(numSamples, 1);
% disp('Generating dataset...');
% for i = 1:numSamples
%     layout = zeros(1,64);
%     stiffIndices = randperm(64, 32);
%     layout(stiffIndices) = 1;
%     X(i,:) = layout;
%     Y(i) = runFEM_layout(layout);
%     if mod(i, 50) == 0, fprintf('Sample %d complete\n', i); end
% end
% save('layout_dataset.mat', 'X', 'Y');