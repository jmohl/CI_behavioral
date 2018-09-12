%% splitting data into train and test sets
%
% -------------------
% Jeff Mohl
% 8/3/18
% -------------------
%
% Description: important to avoid overfitting but splitting the dataset
% into train and test sets, and only comparing the likelihood on test sets.
%
% Methods

% kfold - uses k-fold cross validation where the train set conists of
% (k-1)/k fraction of the data and the test set is the remaining 1/k
% fraction of trials. Trials are selected evenly but not randomly
% throughout dataset. k is specified in options.kfolds and the data is
% returned as an array where each row is a different fold.


%% split data into train and test sets
function [train_inds, test_inds]=get_train_test(data,k)
    test_inds = zeros(k,height(data));
        for i = 1:k
            test_inds(i,i:k:height(data)) = 1;
        end
        train_inds = ~test_inds;
end
