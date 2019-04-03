%% splitting data into train and test sets
%
% -------------------
% Jeff Mohl
% 8/3/18
% updated 3/22/19
% -------------------
%
% Description: important to avoid overfitting but splitting the dataset
% into train and test sets, and only comparing the likelihood on test sets.
%
% Method
% k-fold - uses k-fold cross validation where the train set conists of
% (k-1)/k fraction of the data and the test set is the remaining 1/k
% fraction of trials. k is specified in fitoptions.kfolds and the data is
% returned as 2 cell arrays containing k subsets of tidy_data structures.


%% split data into train and test sets
function [AV_train, AV_test]=get_train_test(data,k_folds)

AV_test = cell(k_folds,1);
AV_train = cell(k_folds,1);

for this_fold = 1:k_folds
    test_inds = logical(zeros(height(data),1));
    test_inds(this_fold:k_folds:height(data)) = 1;
    AV_test{this_fold,1} = data(test_inds,:);
    AV_train{this_fold,1}= data(~test_inds,:);
    
end
