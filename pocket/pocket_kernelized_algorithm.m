prompt = 'SELECT THE DATA SET: ABALONE - SKIN - SHUTTLE: ';
set_name = input(prompt,'s');
switch set_name
        case 'ABALONE'
            data = xlsread('data/abalone.xlsx');
            test_set = data(1:size(data, 1)/2, :); % test set is the first half of the data
            X_test = test_set(:,1:size(test_set, 2)-1); % removes last column(class label)
            y_test = test_set(:,size(test_set, 2):size(test_set, 2)); %last column is the class label
            training_set = data((size(data, 1)/2 + 1):end, :); % training set is the second half of the data
            X = training_set(:,1:size(training_set, 2)-1); % removes last column(class label)
            y = training_set(:,size(training_set, 2):size(training_set, 2)); %last column is the class label
            
            max_step = 10*size(X, 1); 
            
        case 'SKIN'
            data = table2array(readtable('data/skin.txt'));
            test_set = data(1:5000, :); % take the first 5000 as test set
            X_test = test_set(:,1:size(test_set, 2)-1); % removes last column(class label)
            y_test = test_set(:,size(test_set, 2):size(test_set, 2)); %last column is the class label
            y_test = (y_test - 1); % transform 2 in 1 and 1 in 0
            training_set = data(5001:10000, :); % take the second 5000 as training set
            X = training_set(:,1:size(training_set, 2)-1); % removes last column(class label)
            y = training_set(:,size(training_set, 2):size(training_set, 2)); %last column is the class label
            y = (y - 1); % transform 2 in 1 and 1 in 0
            
            max_step = 3*size(X, 1); 
        
        case 'SHUTTLE'
            training_set = table2array(readtable('data/shuttle_training.txt'));
            test_set = table2array(readtable('data/shuttle_test.txt'));
            X_test = test_set(1:5000,1:size(test_set, 2)-1); % take the first 5000 and removes last column (class label)
            y_test = xlsread('data/shuttle_y_test.xlsx');
            X = training_set(1:5000,1:size(training_set, 2)-1); % take the first 5000 and removes last column (class label)
            y = xlsread('data/shuttle_y_training.xlsx');
            
            max_step = 5*size(X, 1); 
        
    otherwise
            fprintf('Error, no such data set!\n');
            return;
end

prompt = 'SELECT A KERNEL: COSINE - POLYNOMIAL2 - YANG: ';
kernel = input(prompt,'s');
if ~(strcmp(kernel,'COSINE') || strcmp(kernel,'POLYNOMIAL2') || strcmp(kernel,'YANG'))    
    fprintf('Error, no such data kernel!\n');
    return;
end
   

X = [ones(size(X, 1), 1) X]; % add column of ones
y = 2*y - 1;  % transform 0 in -1 
X_test = [ones(size(X_test, 1), 1) X_test]; % add column of ones
y_test = 2*y_test - 1;  % transform 0 in -1 

alpha = zeros(size(X,1), 1); % number of times each data point was missclasified

eta = 1;         % the coefficient for the update rule (0 < eta <= 1)
step = 1;
run = 0;
best_run = 0;
max_run = size(X, 1);
alpha_pocket = alpha;


% updates will stop if the number of steps exceeds some maximum number
% or if the run is long enough
i = 1;
while step <= max_step && run < max_run

    if i > size(X,1)  % begin again
        i = 1;
    end
    
    XW = 0;
    for j=1:size(X,1)
        XW = XW + alpha(j) * y(j) * apply_kernel(X(j,:),X(i,:), kernel);
    end
        
    % predict class label for this data point
    y_hat = sign(XW);
    
     if y_hat == 0   % choose 0 as missclassified
        y_hat = -1 * y(i);
    end
    
    if y(i) == y_hat
        run = run + 1;
    else
        if run > best_run
            best_run = run;
            alpha_pocket = alpha;
            run = 0;
        end
        alpha(i) = alpha(i) + 1;
    end
    i = i + 1;
    step = step + 1;
end
if run > best_run
    alpha_pocket = alpha;
end


% Test the weights with our test set
correct_classified = 0;
for i=1:size(X_test,1)
    XW = 0;
    for j=1:size(X,1)
        XW = XW + alpha_pocket(j) * y(j) * apply_kernel(X(j,:),X_test(i,:),kernel);
    end
    
    if sign(XW) == sign(y_test(i))
        correct_classified = correct_classified + 1;
    end
end
classification_accurracy = (correct_classified / size(X_test,1)) * 100;
disp(classification_accurracy);


% Kernel functions
function k = apply_kernel(u, v, type)
    switch type
        case 'POLYNOMIAL2'
            k = compute_polynomial_2_kernel(u, v);
        case 'COSINE'
            k = compute_cosine_kernel(u,v);
        case 'YANG'
            k = compute_yang_kernel(u,v);
        case 'INNER'
            k = compute_inner_product_kernel(u, v);
        otherwise
            fprintf('Error, no such kernel is found!\n')
    end
end

function k = compute_inner_product_kernel(u,v)
    k = u * v';
end

function k = compute_polynomial_2_kernel(u,v)
    k = (u * v'  + 1)^2;
end

function k = compute_cosine_kernel(u,v)
    k = (u * v') / (norm(u) * norm(v));
end

function k = compute_yang_kernel(u,v)
    
    first_sum = 0;
    second_sum = 0;
    for j=1:size(u,2)
        if u(j) >= v(j)
            first_sum = first_sum + (u(j) - v(j));
        else
            second_sum = second_sum + (v(j) - u(j));
        end
    end
    
    sum_maxs = 0;
    for j=1:size(u,2)
        sum_maxs = sum_maxs + max([abs(u(j)), abs(v(j)), abs(u(j) - v(j))]);
    end
    
    k = 1 - sqrt(first_sum^2 + second_sum^2)/sum_maxs;

end
