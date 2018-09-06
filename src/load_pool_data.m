%% load data files and pool across days
%
% -------------------
% Jeff Mohl
% 5/29/18
% -------------------
%
% Description: function loads tidy data files for a given number of days
% (randomly selected, with the option for a seed) and concatenate them into
% a single data table.

function data = load_pool_data(n_pooled_days,seed)

data_dir = dir(['data' '/Juno*.mat']);

if nargin == 2
    rng(seed) %set seed if given
end

rand_days = randsample(length(data_dir),n_pooled_days);

data = [];
for i=1:length(rand_days)
    load(data_dir(rand_days(i)).name); %will be named tidy_data
    data=[data;tidy_data];
end

end