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
function [AV_train, AV_test]=get_train_test(AV_data,method,options)

switch method
    case 'halved'
        for i = 1:length(AV_data)
            train_ind = logical(zeros(length(AV_data{i}),1));
            %randomly pick half the trials in each condition to assign to training
            %dataset, setting those rows to 1
            ntrs = length(AV_data{i});
            train_ind(datasample(1:ntrs,round(ntrs/2),'Replace',false)) = 1;

            AV_train{i} = AV_data{i}(train_ind);
            AV_test{i} = AV_data{i}(~train_ind);
        end
    case 'k-fold'
        k=options.kfolds;
        
        train_ind = logical(zeros(length(AV_data{i}),k));
        
    otherwise
        disp('method not recognized')
end
