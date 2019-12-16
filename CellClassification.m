clear all

load('dataset1.mat');
load('ImageRegistry_v1.mat');
load('LocalizationReference_v1.mat');

%Divide master dataset into training and testing sets
starter = true;
for i = 1:11
    if sum(cellMasterRoster(:, 1) == i) > 0
        subsetRoster = cellMasterRoster(cellMasterRoster(:, 1) == i, :);
        subsetData = cellMasterData(cellMasterRoster(:, 1) == i);
        if starter
            starter = false;
            index = randi(length(subsetRoster(:, 1)));
            testingData = subsetData(index);
            testingRoster = subsetRoster(index, :);
            toKeep = true(length(subsetRoster(:, 1)), 1);
            toKeep(index) = false;
            subsetData = subsetData(toKeep);
            subsetRoster = subsetRoster(toKeep, :);
            for j = 2:round(0.3 * length(subsetRoster(:, 1)))
                index = randi(length(subsetRoster(:, 1)));
                testingData = [testingData; subsetData(index)];
                testingRoster = [testingRoster; subsetRoster(index, :)];
                toKeep = true(length(subsetRoster(:, 1)), 1);
                toKeep(index) = false;
                subsetData = subsetData(toKeep);
                subsetRoster = subsetRoster(toKeep, :);
            end
            trainingData = subsetData;
            trainingRoster = subsetRoster;
        else
            for j = 1:round(0.3 * length(subsetRoster(:, 1)))
                index = randi(length(subsetRoster(:, 1)));
                testingData = [testingData; subsetData(index)];
                testingRoster = [testingRoster; subsetRoster(index, :)];
                toKeep = true(length(subsetRoster(:, 1)), 1);
                toKeep(index) = false;
                subsetData = subsetData(toKeep);
                subsetRoster = subsetRoster(toKeep, :);
            end
            trainingData = [trainingData; subsetData];
            trainingRoster = [trainingRoster; subsetRoster];
        end
    end
end

save('dataset1_breakdown.mat', 'trainingData', 'trainingRoster', 'testingData', 'testingRoster');



%LBP and HOG Feature Extraction
load('dataset1_breakdown.mat');

trainingLBPFeatures = zeros(length(trainingData), 10);
trainingHOGFeatures = zeros(length(trainingData), 324);
for i = 1:length(trainingData)
    % test = trainingData(1).Image;
    % imshow(test)
    test2 = uint8(zeros(512, 672));
    test2(trainingData(i).PixelIdxList) = trainingData(i).PixelValues;
    % [width, height] = size(trainingData(1).Image);
    % test2 = imcrop(test2, [min(trainingData(1).PixelList(:, 1)) min(trainingData(1).PixelList(:, 2)) width height]);
    test2 = test2(min(trainingData(i).PixelList(:, 2)):max(trainingData(i).PixelList(:, 2)), min(trainingData(i).PixelList(:, 1)):max(trainingData(i).PixelList(:, 1)));
    wdiff = 35 - length(test2(1, :));
    hdiff = 35 - length(test2(:, 1));
    if wdiff > 0
        leftadd = round(wdiff / 2);
        rightadd = wdiff - leftadd;
        test2 = [zeros(length(test2(:, 1)), leftadd), test2, zeros(length(test2(:, 1)), rightadd)];
    end
    if hdiff > 0
        topadd = round(hdiff / 2);
        bottomadd = hdiff - topadd;
        test2 = [zeros(topadd, length(test2(1, :))); test2; zeros(bottomadd, length(test2(1, :)))];
    end

    trainingLBPFeatures(i, :) = extractLBPFeatures(test2, 'Upright', false);
    trainingHOGFeatures(i, :) = extractHOGFeatures(test2);
    
    if mod(i, 100) == 0
        disp(i)
    end
end

testingLBPFeatures = zeros(length(testingData), 10);
testingHOGFeatures = zeros(length(testingData), 324);
for i = 1:length(testingData)
    % test = trainingData(1).Image;
    % imshow(test)
    test2 = uint8(zeros(512, 672));
    test2(testingData(i).PixelIdxList) = testingData(i).PixelValues;
    % [width, height] = size(trainingData(1).Image);
    % test2 = imcrop(test2, [min(trainingData(1).PixelList(:, 1)) min(trainingData(1).PixelList(:, 2)) width height]);
    test2 = test2(min(testingData(i).PixelList(:, 2)):max(testingData(i).PixelList(:, 2)), min(testingData(i).PixelList(:, 1)):max(testingData(i).PixelList(:, 1)));
    wdiff = 35 - length(test2(1, :));
    hdiff = 35 - length(test2(:, 1));
    if wdiff > 0
        leftadd = round(wdiff / 2);
        rightadd = wdiff - leftadd;
        test2 = [zeros(length(test2(:, 1)), leftadd), test2, zeros(length(test2(:, 1)), rightadd)];
    end
    if hdiff > 0
        topadd = round(hdiff / 2);
        bottomadd = hdiff - topadd;
        test2 = [zeros(topadd, length(test2(1, :))); test2; zeros(bottomadd, length(test2(1, :)))];
    end

    testingLBPFeatures(i, :) = extractLBPFeatures(test2, 'Upright', false);
    testingHOGFeatures(i, :) = extractHOGFeatures(test2);
    
    if mod(i, 100) == 0
        disp(i)
    end
end


%Run KNN for LBP Features
locRatio = [1, 1, 10, 2, 3, 1, 1, 6, 3, 1, 1];

knn_training_data = zeros(size(trainingLBPFeatures));
for i = 1:10
    knn_training_data(:, i) = trainingLBPFeatures(:, i) - min(trainingLBPFeatures(:, i));
    knn_training_data(:, i) = knn_training_data(:, i) / max(knn_training_data(:, i));
end
knn_training_labels = categorical(trainingRoster(:, 1));
knn_testing_data = zeros(size(testingLBPFeatures));
for i = 1:10
    knn_testing_data(:, i) = testingLBPFeatures(:, i) - min(testingLBPFeatures(:, i));
    knn_testing_data(:, i) = knn_testing_data(:, i) / max(knn_testing_data(:, i));
end
knn_testing_labels = categorical(testingRoster(:, 1));

%Scanning the number of neighbors hyperparameter to identify best value to
%use
knn_predictions = zeros(length(knn_testing_labels), 100);
numneigh = zeros(1, 100);
overall_accuracy = zeros(1, 100);
knn_testing_labels_temp = double(knn_testing_labels);
for i = 1:100
    knn_model = fitcknn(knn_training_data, knn_training_labels, 'Prior', locRatio./sum(locRatio), 'Numneighbors', i);

    knn_predictions(:, i) = predict(knn_model, knn_testing_data);

    cmatrix = confusionmat(knn_testing_labels_temp, knn_predictions(:, i));
%     confusionchart(knn_testing_labels, knn_predictions);

    good = 0;
    for j = 1:11
        good = good + cmatrix(j, j);
    end
    numneigh(i) = i;
    overall_accuracy(i) = good / sum(sum(cmatrix));
end

figure();
plot(numneigh, overall_accuracy, 'LineWidth', 2);
xlabel('Number of neighbors used in KNN Classifier');
ylabel('Fraction of cell labels correctly predicted');

figure();
confusionchart(knn_testing_labels_temp, knn_predictions(:, 64));

%Using matlab's built-in optimization of hyperparameters
knn_model2 = fitcknn(knn_training_data, knn_training_labels, 'Prior', locRatio./sum(locRatio), 'OptimizeHyperparameters', 'auto', 'HyperparameterOptimizationOptions', struct('AcquisitionFunctionName','expected-improvement-plus'));
knn_predictions2 = predict(knn_model2, knn_testing_data);
confusionmat(knn_testing_labels, knn_predictions2)
figure();
confusionchart(knn_testing_labels, knn_predictions2);


%Compare both sets of predictions to true values
cmatrix1 = confusionmat(knn_testing_labels_temp, knn_predictions(:, 64));
cmatrix2 = confusionmat(knn_testing_labels, knn_predictions2);
good1 = zeros(1, 11);
good2 = zeros(1, 11);
for j = 1:11
    good1(j) = cmatrix1(j, j);
    good2(j) = cmatrix2(j, j);
end
figure();
bar([1 2 3 4 5 6 7 8 9 10 11], [sum(cmatrix1, 2), good1', good2']);
xlabel('Localization class');
ylabel('Cell Number');
legend('Actual cells', 'Correct predictions - non-optimized', 'Correct predictions - optimized');






%Rational binary decision tree design
rational_training_data(1).PixelValues = double(trainingData(1).PixelValues)./double(max(trainingData(1).PixelValues));
for i = 2:length(trainingData)
    rational_training_data(i).PixelValues = double(trainingData(i).PixelValues)./double(max(trainingData(i).PixelValues));
end
rational_training_labels = trainingRoster(:, 1);

%Metric 1
%Nucleolus and nucleus
dataOfInterest = rational_training_data(rational_training_labels == 8 | rational_training_labels == 7);
labelsOfInterest = rational_training_labels(rational_training_labels == 8 | rational_training_labels == 7);
scores = zeros(100);
intensities = zeros(100, 1);
areas = zeros(100, 1);
for i = 0.01:0.01:1
    intensities(int16(i * 100)) = i;
    fractions = zeros(length(dataOfInterest), 1);
    for j = 1:length(dataOfInterest)
        fractions(j) = sum(dataOfInterest(j).PixelValues >= i) / length(dataOfInterest(j).PixelValues);
    end
    for j = 0.01:0.01:1
        areas(int16(j * 100)) = j;
        state = fractions >= j;
        output = zeros(length(fractions), 1);
        for k = 1:length(state)
            if state(k)
                output(k) = 8;
            else
                output(k) = 7;
            end
        end
        outcome = confusionmat(labelsOfInterest, output);
        %below, consider nuclear localization to be "positive"
        TP = outcome(2, 2);
        TN = outcome(1, 1);
        FP = outcome(1, 2);
        FN = outcome(2, 1);
        scores(int16(i * 100), int16(j * 100)) = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN));
    end
end

figure();
mesh(areas, intensities, scores);
[temp, I_i] = max(scores);
[~, I_j] = max(temp);
threshold1 = [I_i(I_j) / 100, I_j / 100];

%Punctate and mitochondria
dataOfInterest = rational_training_data(rational_training_labels == 5 | rational_training_labels == 9);
labelsOfInterest = rational_training_labels(rational_training_labels == 5 | rational_training_labels == 9);
scores = zeros(100);
intensities = zeros(100, 1);
areas = zeros(100, 1);
for i = 0.01:0.01:1
    intensities(int16(i * 100)) = i;
    fractions = zeros(length(dataOfInterest), 1);
    for j = 1:length(dataOfInterest)
        fractions(j) = sum(dataOfInterest(j).PixelValues >= i) / length(dataOfInterest(j).PixelValues);
    end
    for j = 0.01:0.01:1
        areas(int16(j * 100)) = j;
        state = fractions >= j;
        output = zeros(length(fractions), 1);
        for k = 1:length(state)
            if state(k)
                output(k) = 9;
            else
                output(k) = 5;
            end
        end
        outcome = confusionmat(labelsOfInterest, output);
        %below, consider punctate localization to be "positive"
        TP = outcome(2, 2);
        TN = outcome(1, 1);
        FP = outcome(1, 2);
        FN = outcome(2, 1);
        scores(int16(i * 100), int16(j * 100)) = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN));
    end
end

figure();
mesh(areas, intensities, scores);
[temp, I_i] = max(scores);
[~, I_j] = max(temp);
threshold2 = [I_i(I_j) / 100, I_j / 100];

%Nuclear periphery, nucleolus, nucleus, and budneck from vacuole and
%vacuole membrane
dataOfInterest = rational_training_data(rational_training_labels == 1 | rational_training_labels == 6 | rational_training_labels == 7 | rational_training_labels == 8 | rational_training_labels == 10 | rational_training_labels == 11);
labelsOfInterest = rational_training_labels(rational_training_labels == 1 | rational_training_labels == 6 | rational_training_labels == 7 | rational_training_labels == 8 | rational_training_labels == 10 | rational_training_labels == 11);
for i = 1:length(labelsOfInterest)
    if labelsOfInterest(i) == 1 || labelsOfInterest(i) == 6 || labelsOfInterest(i) == 7 || labelsOfInterest(i) == 8
        labelsOfInterest(i) = 1;
    else
        labelsOfInterest(i) = 2;
    end
end
scores = zeros(100);
intensities = zeros(100, 1);
areas = zeros(100, 1);
for i = 0.01:0.01:1
    intensities(int16(i * 100)) = i;
    fractions = zeros(length(dataOfInterest), 1);
    for j = 1:length(dataOfInterest)
        fractions(j) = sum(dataOfInterest(j).PixelValues >= i) / length(dataOfInterest(j).PixelValues);
    end
    for j = 0.01:0.01:1
        areas(int16(j * 100)) = j;
        state = fractions >= j;
        output = zeros(length(fractions), 1);
        for k = 1:length(state)
            if state(k)
                output(k) = 2;
            else
                output(k) = 1;
            end
        end
        outcome = confusionmat(labelsOfInterest, output);
        %below, consider punctate localization to be "positive"
        TP = outcome(2, 2);
        TN = outcome(1, 1);
        FP = outcome(1, 2);
        FN = outcome(2, 1);
        scores(int16(i * 100), int16(j * 100)) = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN));
    end
end

figure();
mesh(areas, intensities, scores);
[temp, I_i] = max(scores);
[~, I_j] = max(temp);
threshold3 = [I_i(I_j) / 100, I_j / 100];


%Metric 2
%All classes
labelsOfInterest = trainingRoster(:, 1);
for i = 1:length(labelsOfInterest)
    if labelsOfInterest(i) == 2 || labelsOfInterest(i) == 3 || labelsOfInterest(i) == 4
        labelsOfInterest(i) = 1;
    else
        labelsOfInterest(i) = 2;
    end
end
radii = sqrt([trainingData.Area]./pi);
centroids = zeros(length(trainingData), 2);
weighted = zeros(length(trainingData), 2);
for i = 1:length(trainingData)
    centroids(i, :) = trainingData(i).Centroid;
    weighted(i, :) = trainingData(i).WeightedCentroid;
end
distances = sqrt(((centroids(:, 1) - weighted(:, 1)).^2) + ((centroids(:, 2) - weighted(:, 2)).^2))./radii';
scores = zeros(100, 1);
limits = linspace(min(distances), max(distances));
for i = 1:100
    state = distances >= limits(i);
    output = zeros(length(distances), 1);
    for j = 1:length(state)
        if state(j)
            output(j) = 2;
        else
            output(j) = 1;
        end
    end
    outcome = confusionmat(labelsOfInterest, output);
    %below, consider shifted localization to be "positive"
    TP = outcome(2, 2);
    TN = outcome(1, 1);
    FP = outcome(1, 2);
    FN = outcome(2, 1);
    scores(i) = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN));
end

figure();
plot(limits, scores);
[~, I] = max(scores);
threshold4 = I;
