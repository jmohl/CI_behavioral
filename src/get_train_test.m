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
% halved - equally splits each condition type in half, and return the
% resulting data structures as 2 1x20 cell arrays, matching the input
% structure
% k-fold - uses k-fold cross validation where the train set conists of
% (k-1)/k fraction of the data and the test set is the remaining 1/k
% fraction of triaks. k is specified in options.kfolds and the data is
% returned as 2 kx20 cell arrays.


%% split data into train and test sets
function [AV_train, AV_test]=get_train_test(data,k_folds, AV_pairs)

AV_test = cell(k_folds, length(AV_pairs));
AV_train = cell(k_folds,length(AV_pairs));

for this_fold = 1:k_folds
    test_inds = logical(zeros(height(data),1));
    test_inds(this_fold:k_folds:height(data)) = 1;
    test_data = data(test_inds,:);
    train_data = data(~test_inds,:);
    %split into cell array for faster fitting code
    for this_pair = 1:length(AV_pairs)
        A_tar = AV_pairs(this_pair,1);
        V_tar = AV_pairs(this_pair,2);
        try %keeps from failing when there is no
        test_endpoints = vertcat(test_data(test_data.A_tar == A_tar & test_data.V_tar == V_tar,:).valid_endpoints{:,1});
        test_endpoints = test_endpoints(:,1); %only x component
        AV_test{this_fold,this_pair} = test_endpoints;
        end
        try
        train_endpoints = vertcat(train_data(train_data.A_tar == A_tar & train_data.V_tar == V_tar,:).valid_endpoints{:,1});
        train_endpoints = train_endpoints(:,1); %only x component
        AV_train{this_fold,this_pair} = train_endpoints;
        end
    end
    
end
