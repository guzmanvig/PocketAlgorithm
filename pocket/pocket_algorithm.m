% data = xlsread('debt.xls');
% % Test set as the first half of data
% test_set = data(1:floor(size(data, 1)/2), :);
% X_test = test_set(:,2:size(test_set, 2)-1); % removes first column (id) and last (class label)
% y_test = test_set(:,size(test_set, 2):size(test_set, 2)); %last column is the class label
% % Training set as the second half of data
% training_set = data((floor(size(data, 1)/2) + 1):end, :);
% X = training_set(:,2:size(training_set, 2)-1); % removes first column (id) and last (class label)
% y = training_set(:,size(training_set, 2):size(training_set, 2)); %last column is the class label

% data = xlsread('iris_shuffled.xlsx');
% % Test set as the first half of the shuffled data
% test_set = data(1:size(data, 1)/2, :);
% X_test = test_set(:,1:size(test_set, 2)-1); % removes last column(class label)
% y_test = test_set(:,size(test_set, 2):size(test_set, 2)); %last column is the class label
% % Training set as the second half of the shuffled data
% training_set = data((size(data, 1)/2 + 1):end, :);
% X = training_set(:,1:size(training_set, 2)-1); % removes last column(class label)
% y = training_set(:,size(training_set, 2):size(training_set, 2)); %last column is the class label

% data = table2array(readtable('breast_cancer_no_nan.txt'));
% % Test set as the first half of data
% test_set = data(1:floor(size(data, 1)/2), :);
% X_test = test_set(:,2:size(test_set, 2)-1); % removes first column (id) and last (class label)
% y_test = test_set(:,size(test_set, 2):size(test_set, 2)); %last column is the class label
% y_test = (y_test - 2) / 2; % transform 2 in 0 and 4 in 1
% % Training set as the second half of data
% training_set = data((floor(size(data, 1)/2) + 1):end, :);
% X = training_set(:,2:size(training_set, 2)-1); % removes first column (id) and last (class label)
% y = training_set(:,size(training_set, 2):size(training_set, 2)); %last column is the class label
% y = (y - 2) / 2; % transform 2 in 0 and 4 in 1

% training_set = table2array(readtable('shuttle_training.txt'));
% test_set = table2array(readtable('shuttle_test.txt'));
% X_test = test_set(1:5000,1:size(test_set, 2)-1); % removes last column (class label)
% y_test = xlsread('shuttle_y_test.xlsx');
% X = training_set(1:5000,1:size(training_set, 2)-1); % removes first column (id) and last (class label)
% y = xlsread('shuttle_y_training.xlsx');

% data = xlsread('wholesale_data.csv');
% % Test set as the first half of data
% test_set = data(1:floor(size(data, 1)/2), :);
% X_test = test_set(:,2:size(test_set, 2)); % removes first column(class label)
% y_test = test_set(:,1); %first column is the class label
% y_test = (y_test - 1); % transform 2 in 1 and 1 in 0
% % Training set as the second half of data
% training_set = data((floor(size(data, 1)/2) + 1):end, :);
% X = training_set(:,2:size(test_set, 2)); % removes first column(class label)
% y = training_set(:,1); %first column is the class label
% y = (y - 1); % transform 2 in 1 and 1 in 0

% data = xlsread('HTRU.csv');
% % Test set as the first half of the shuffled data
% test_set = data(1:size(data, 1)/2, :);
% X_test = test_set(:,1:size(test_set, 2)-1); % removes last column(class label)
% y_test = test_set(:,size(test_set, 2):size(test_set, 2)); %last column is the class label
% % Training set as the second half of the shuffled data
% training_set = data((size(data, 1)/2 + 1):end, :);
% X = training_set(:,1:size(training_set, 2)-1); % removes last column(class label)
% y = training_set(:,size(training_set, 2):size(training_set, 2)); %last column is the class label

data = table2array(readtable('Skin_cut.txt'));
% Test set as the first half of the shuffled data
%test_set = data(1:size(data, 1)/2, :);
test_set = data(1:5000, :);
X_test = test_set(:,1:size(test_set, 2)-1); % removes last column(class label)
y_test = test_set(:,size(test_set, 2):size(test_set, 2)); %last column is the class label
y_test = (y_test - 1); % transform 2 in 1 and 1 in 0
% Training set as the second half of the shuffled data
%training_set = data((size(data, 1)/2 + 1):end, :);
training_set = data(5001:10000, :);
X = training_set(:,1:size(training_set, 2)-1); % removes last column(class label)
y = training_set(:,size(training_set, 2):size(training_set, 2)); %last column is the class label
y = (y - 1); % transform 2 in 1 and 1 in 0

% data = xlsread('position_modified_cut.xlsx');
% % Test set as the first half of the shuffled data
% %test_set = data(1:size(data, 1)/2, :);
% test_set = data(1:5000, :);
% X_test = test_set(:,1:size(test_set, 2)-1); % removes last column(class label)
% y_test = test_set(:,size(test_set, 2):size(test_set, 2)); %last column is the class label
% % Training set as the second half of the shuffled data
% %training_set = data((size(data, 1)/2 + 1):end, :);
% training_set = data(5001:10000, :);
% X = training_set(:,1:size(training_set, 2)-1); % removes last column(class label)
% y = training_set(:,size(training_set, 2):size(training_set, 2)); %last column is the class label

% data = xlsread('abalone_01_shuffled.xlsx');
% %suffled_data = data(randperm(size(data, 1)), :);
% % Test set as the first half of the shuffled data
% test_set = data(1:size(data, 1)/2, :);
% X_test = test_set(:,1:size(test_set, 2)-1); % removes last column(class label)
% y_test = test_set(:,size(test_set, 2):size(test_set, 2)); %last column is the class label
% % Training set as the second half of the shuffled data
% training_set = data((size(data, 1)/2 + 1):end, :);
% X = training_set(:,1:size(training_set, 2)-1); % removes last column(class label)
% y = training_set(:,size(training_set, 2):size(training_set, 2)); %last column is the class label

% % Test set as the first half of the shuffled data
% test_set = xlsread('contraceptive_test.xlsx');
% training_set = xlsread('contraceptive_train.xlsx');
% X_test = test_set(:,1:size(test_set, 2)-1); % removes last column(class label)
% y_test = test_set(:,size(test_set, 2):size(test_set, 2)); %last column is the class label
% % Training set as the second half of the shuffled data
% X = training_set(:,1:size(training_set, 2)-1); % removes last column(class label)
% y = training_set(:,size(training_set, 2):size(training_set, 2)); %last column is the class label


X = [ones(size(X, 1), 1) X]; % add column of ones
y = 2*y - 1;  % transform 0 in -1 
X_test = [ones(size(X_test, 1), 1) X_test]; % add column of ones
y_test = 2*y_test - 1;  % transform 0 in -1 

w = zeros(size(X,2), 1); %initialize weights as a column of 0 of d+1 dimension

% initialize stuff
max_step = 5*size(X, 1); % maximum number of iterations
eta = 1;         % the coefficient for the update rule (0 < eta <= 1)
step = 1;
run = 0;
best_run = 0;
max_run = size(X, 1);
w_pocket = w;

% updates will stop if the number of steps exceeds some maximum number
% or if the run is long enough
i = 1;
while step <= max_step && run < max_run

    if i > size(X,1)  % begin again
        i = 1;
    end
    % predict class label for this data point
    y_hat = sign(X(i, :) * w);
    
    if y_hat == 0   % choose 0 as missclassified
        y_hat = -1 * y(i);
    end
    
    if y(i) == y_hat
        run = run + 1;
    else
        if run > best_run
            best_run = run;
            w_pocket = w;
            run = 0;
        end
        w = w + eta * 0.5 * (y(i) - y_hat) * X(i, :)';
    end
  
    i = i + 1;
    step = step + 1;
end
if run > best_run
    w_pocket = w;
end

% Test the weights with our test set
correct_classified = 0;
for i=1:size(X_test,1)
    if sign(X_test(i,:) * w_pocket) == sign(y_test(i))
        correct_classified = correct_classified + 1;
    end
end
classification_accurracy = (correct_classified / size(X_test,1)) * 100;
disp(classification_accurracy);
