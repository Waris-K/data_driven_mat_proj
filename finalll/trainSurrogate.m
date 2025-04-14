function [net, testMSE, R2] = trainSurrogate()
% Trains a neural network on the dataset and evaluates its performance

load('fem_dataset.mat', 'data');
X = data(:, 1:end-1);
Y = data(:, end);

%% Split data
rng("default") % For reproducibility
c = cvpartition(length(Y), "Holdout", 0.20);
trainingIdx = training(c);
XTrain = X(trainingIdx,:);
YTrain = Y(trainingIdx);
testIdx = test(c);
XTest = X(testIdx,:);
YTest = Y(testIdx);

%% Train the surrogate neural network
net = fitrnet(XTrain, YTrain, ...
    'LayerSizes', [30 10], ...
    'Activations', 'relu', ...
    'Standardize', true);

save('trained_surrogate.mat', 'net');

%% Test MSE and predictions
testMSE = loss(net, XTest, YTest);
testPredictions = predict(net, XTest);

%% Plot True vs Predicted
figure;
scatter(YTest, testPredictions, 'filled')
hold on
plot(YTest, YTest, 'r--', 'LineWidth', 2)
hold off
xlabel("True SE")
ylabel("Predicted SE")
title("True vs Predicted SE")
grid on

%% R^2 calculation
SS_res = sum((YTest - testPredictions).^2);
SS_tot = sum((YTest - mean(YTest)).^2);
R2 = 1 - SS_res / SS_tot;
fprintf("Test MSE = %.4f\n", testMSE);
fprintf("R^2 on test set = %.4f\n", R2);

% %% Visualize network
% figure;
% plot(net)
% title('Neural Network Architecture')

end
