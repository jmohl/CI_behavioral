%% master script for running Causal inference analysis and figure generation
%
% -------------------
% Jeff Mohl
% 5/29/18
% -------------------
%
% Description: this script will allow running of each step of the
% behavioral CI project, from loading data, fitting model, and generating
% figures


%% hardcoded parameters
local_directory = 'C:\Users\jtm47\Documents\Projects\CI_behavioral\';
n_pooled_days = 15;
seed = 'default';

%adding paths
cd(local_directory)
addpath('src','results','data')


%% load data
data = load_pool_data(n_pooled_days,seed); %data is pooled across N randomly selected days, yielding a single tidy data table 

[AV_endpoint_data,AV_tar_pairs] = get_endpoint_array(data,'AV',1,100); %get only relevant saccade endpoints, for AV trials, used for fitting model

% [V_endpoint_data,V_tars] = get_endpoint_array(data,'V',1,100); %get only relevant saccade endpoints, for AV trials
% [A_endpoint_data,A_tars] = get_endpoint_array(data,'A',1,100); %get only relevant saccade endpoints, for AV trials

%% fit CI model
% load initial parameters
load('ini_params.mat')
options = optimset('MaxFunEvals',1000,'MaxIter',2000,'PlotFcns',@optimplotfval);
free_params_vector = cell2mat(struct2cell(free_params)); %need to cast params as vector for fminsearch

% run fminsearch to fit parameters
CI_minsearch = @(free_params_vector)get_nll_CI(AV_tar_pairs,AV_endpoint_data,fixed_params,free_params_vector);
[fit_params_vector,min_nll_CI,~,~] = fminsearch(CI_minsearch,free_params_vector,options);

%put params back in structure format
fit_params_CI = free_params;
names = fieldnames(free_params);
for i=1:length(names)
fit_params_CI.(names{i})= fit_params_vector(i);
end

%% fit alternative models

%% model comparison

%% generate figures
%plotting demos
A_tar = -24;
V_tar = -6;
plot_fit_dists(data,A_tar,V_tar,fit_params_CI,fixed_params,-40:40);

