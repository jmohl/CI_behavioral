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

% load initial parameters
load('ini_params.mat')


%% load data
data = load_pool_data(n_pooled_days,seed); %data is pooled across N randomly selected days, yielding a single tidy data table 
%data = load('H02_AVD2_2018_07_24_tidy.mat'); 
%data=data.this_tidy;
[AV_endpoint_data,AV_tar_pairs] = get_endpoint_array(data,'AV',1,100); %get only relevant saccade endpoints, for AV trials, used for fitting model
[V_endpoint_data,fixed_params.V_tars] = get_endpoint_array(data,'V',1,100); %get only relevant saccade endpoints, for V trials
[A_endpoint_data,fixed_params.A_tars] = get_endpoint_array(data,'A',1,100); %get only relevant saccade endpoints, for A trials

%% get single modality estimates
[fixed_params.V_mu,fixed_params.V_sig] = get_unimodal_est(V_endpoint_data);
[fixed_params.A_mu,fixed_params.A_sig] = get_unimodal_est(A_endpoint_data);


%% fit CI model
%sandboxing stuff, for temporary use
%restricted_trials = AV_tar_pairs(:,1) ~= 6 & AV_tar_pairs(:,1) ~= -6;
%
options = optimset('MaxFunEvals',10000,'MaxIter',20000,'PlotFcns',@optimplotfval);
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

%% Alternative models

% fully segregated model

% fully integrated model

% probability matching instead of posterior reweighted



%CI model with target locs and sigmas determined by unimodal fits
options = optimset('MaxFunEvals',1000,'MaxIter',2000,'PlotFcns',@optimplotfval);
free_params_vector = cell2mat(struct2cell(free_params)); %need to cast params as vector for fminsearch
free_params_vector = free_params_vector(4:end);

% run fminsearch to fit parameters
CI_fAV_minsearch = @(free_params_vector)get_nll_CI_fixedAV(AV_tar_pairs,AV_endpoint_data,fixed_params,free_params_vector);
[fit_params_vector_fAV,min_nll_CI_fAV,~,~] = fminsearch(CI_fAV_minsearch,free_params_vector,options);

%model does not converge, ends up with extremely large prior sig and
%extremely low p_common indicating that these variables are adding very
%little. Might make sense to introduce a bias term or something.

%put params back in structure format
fit_params_CI_fAV = free_params;
names = fieldnames(free_params);
for i=1:length(names)
fit_params_CI_fAV.(names{i})= fit_params_vector(i);
end


%% model comparison

%% generate figures
%plotting demos
plot_dir = 'results\fit_dists';
for i = 1:size(AV_tar_pairs,1)
A_tar = AV_tar_pairs(i,1);
V_tar =  AV_tar_pairs(i,2);
plot_fit_dists(data,A_tar,V_tar,fit_params_CI,fixed_params,-40:40);
saveas(gcf,sprintf('%s\\%dA%dV',plot_dir,A_tar,V_tar),'png')
end
