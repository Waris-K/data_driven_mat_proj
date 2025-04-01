clear; clc; close all; 

%% Load the inputs and labels
load('displacement.mat')                % displacements are the lables (outputs)
load('material_distributionC.mat')      % material distribution is the data (input)
%numNodes = 9*9;
%% Splitting testing and training data

displacement = displacement(1:2:end,:); 
outputSize = 81; 

numSamples = size(material_distribution, 3);  % Total number of samples (e.g., 1000)
numTrain = floor(0.8 * numSamples);  % 80% for training
numValidation = numSamples - numTrain;  % 20% for validation

% Split the data into training and validation sets
trainData = material_distribution(:, :, 1:numTrain);  % Training data
trainLabels = displacement(:, 1:numTrain);  % Training labels (displacements)

valData = material_distribution(:, :, numTrain+1:end);  % Validation data
valLabels = displacement(:, numTrain+1:end);  % Validation labels


% reshaping the data and labels to match the expected inputs/outputs
trainData = reshape(trainData, [8, 8, 1, numTrain]);
valData = reshape(valData, [8, 8, 1, numValidation]);

% Reshape displacement labels to match network output format
trainLabels = reshape(trainLabels, [1, 1, outputSize, numTrain]);
valLabels = reshape(valLabels, [1, 1, outputSize, numValidation]);

%% CNN Architecture

% % Define CNN architecture
% layers = [
%     imageInputLayer([8 8 1])  % Input layer for 8x8 material distribution matrix
%     convolution2dLayer(3, 16, 'Padding', 'same')  % Convolutional layer with 16 filters of size 3x3
%     reluLayer  % ReLU activation function
%     convolution2dLayer(3, 32, 'Padding', 'same')  % Convolutional layer with 32 filters of size 3x3
%     reluLayer  % ReLU activation function
%     fullyConnectedLayer(162)  % Fully connected layer with 162 neurons matches output labels
% %     reluLayer  % ReLU activation function
% %     fullyConnectedLayer(2 * numNodes)  % Output layer (displacement values for Ux and Uy)
%     regressionLayer];  % Regression output

layers = [
    imageInputLayer([8 8 1])  
    convolution2dLayer(3, 16, 'Padding', 'same')  
    reluLayer  
    maxPooling2dLayer(2, 'Stride', 2)  % Pooling to downsample
    convolution2dLayer(3, 32, 'Padding', 'same')  
    reluLayer  
    maxPooling2dLayer(2, 'Stride', 2)  
    fullyConnectedLayer(512)  
    reluLayer
    fullyConnectedLayer(outputSize)  
    regressionLayer];  


% Set the options for training with validation data
options = trainingOptions('adam', ...
    'MaxEpochs', 100, ...
    'MiniBatchSize', 32, ...
    'InitialLearnRate', 1e-3, ...
    'Shuffle', 'every-epoch', ...
    'ValidationData', {valData, valLabels},...  % Provide validation data
    'ValidationFrequency', 30,...  % Frequency of validation during training
    'ValidationPatience', 5,...  % Number of epochs to wait for improvement before stopping
    'Plots', 'training-progress', ...
    'Verbose', false);  % Disable verbose output for training

%% Train the network
[netX,info] = trainNetwork(trainData, trainLabels, layers, options);


%% saving the training

% Save the trained network
save('trainedCNNModelxOnly.mat', 'netX');

%% load the network
load('trainedCNNModelxOnly.mat')

%% Evaluate the Validation

% Evaluate the model on the validation data
predictedValidationLabels = predict(netX, valData);

% Compare the predicted values with the actual labels (displacements)
validationError = mean(abs(predictedValidationLabels - valLabels), 'all');  % Mean absolute error

disp(['Validation Error: ', num2str(validationError)]);


%% actual vs predicted
actualDisplacements = squeeze(valLabels);  % Removes singleton dimensions
actualDisplacements = actualDisplacements'; 

predictedDisplacements = predictedValidationLabels;

figure;
scatter(actualDisplacements, predictedDisplacements, 'filled');
hold on;
plot(linspace(min(actualDisplacements(:)), max(actualDisplacements(:)), 100), ...
     linspace(min(actualDisplacements(:)), max(actualDisplacements(:)), 100), ...
     'r--', 'LineWidth', 2); % y=x line (red dashed)
xlabel('Actual Displacement');
ylabel('Predicted Displacement');
title('Actual vs. Predicted Displacement');
grid on;
axis equal; % Ensure 1:1 scaling


%% visualizing the loss in training and validation

% % You can access the training progress in the trainingHistory variable
% history = net.trainingHistory;
% 
% % Plot training and validation losses
% figure;
% plot(history.TrainingLoss, 'b', 'LineWidth', 2);  % Training loss (blue)
% hold on;
% plot(history.ValidationLoss, 'r', 'LineWidth', 2);  % Validation loss (red)
% legend('Training Loss', 'Validation Loss');
% xlabel('Epoch');
% ylabel('Loss');
% title('Training and Validation Loss');
% grid on;
