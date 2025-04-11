
function trainSurrogate()
% Trains a neural network on the dataset

load('fem_dataset.mat', 'data');
X = data(:, 1:end-1);
Y = data(:, end);

%%
rng("default") % For reproducibility of the partition
c = cvpartition(length(Y),"Holdout",0.20);
trainingIdx = training(c); % Indices for the training set
XTrain = X(trainingIdx,:);
YTrain = Y(trainingIdx);
testIdx = test(c); % Indices for the test set
XTest = X(testIdx,:);
YTest = Y(testIdx);
%%

net = fitrnet(XTrain, YTrain, 'LayerSizes', [30 10], 'Activations','relu', 'Standardize', true);
save('trained_surrogate.mat', 'net');

testMSE = loss(net,XTest,YTest)

testPredictions = predict(net,XTest);
plot(YTest,testPredictions,".")
hold on
plot(YTest,YTest)
hold off
xlabel("True se")
ylabel("Predicted se")
end
