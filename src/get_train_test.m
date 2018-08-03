%% splitting data into train and test sets
%
% -------------------
% Jeff Mohl
% 8/3/18
% -------------------
%
% Description: important to avoid overfitting but splitting the dataset
% into train and test sets, and only comparing the likelihood on test sets.
% Currently I am doing this by splitting the dataset in half randomly. In
% the future this could be improved by making several train/test splits and
% using some kind of folded cross validation.


%% split data into train and test sets
function [AV_train, AV_test]=get_train_test(AV_data)
for i = 1:length(AV_data)
    train_ind = logical(zeros(length(AV_data{i}),1));
    %randomly pick half the trials in each condition to assign to training
    %dataset, setting those rows to 1
    ntrs = length(AV_data{i});
    train_ind(datasample(1:ntrs,round(ntrs/2),'Replace',false)) = 1;
    
    AV_train{i} = AV_data{i}(train_ind);
    AV_test{i} = AV_data{i}(~train_ind);
end
end
