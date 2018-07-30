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

%Todo:
% - correct for multiplicative and additive biases, based on visual (accounts
%for improperly calibrated eye tracker)
% - figure out if doing anything about sigmoidal aud loc bias


%% hardcoded parameters
local_directory = 'C:\Users\jtm47\Documents\Projects\CI_behavioral\';
n_pooled_days = 15;
seed = 'default';
make_plots = 0;

%adding paths
cd(local_directory)
addpath('src','results','data')

% load initial parameters
load('ini_params.mat')


%% load data
%data = load_pool_data(n_pooled_days,seed); %data is pooled across N randomly selected days, yielding a single tidy data table 
data = load('H02_AVD2_2018_07_24_tidy.mat'); %running on a single human data
data=data.this_tidy;
[AV_endpoint_data,AV_tar_pairs] = get_endpoint_array(data,'AV',1,100); %get only relevant saccade endpoints, for AV trials, used for fitting model
[V_endpoint_data,fixed_params.V_tars] = get_endpoint_array(data,'V',1,100); %get only relevant saccade endpoints, for V trials
[A_endpoint_data,fixed_params.A_tars] = get_endpoint_array(data,'A',1,100); %get only relevant saccade endpoints, for A trials

%% get single modality estimates
[fixed_params.V_mu,fixed_params.V_sig] = get_unimodal_est(V_endpoint_data);
[fixed_params.A_mu,fixed_params.A_sig] = get_unimodal_est(A_endpoint_data);


%% split data into train and test sets
for i = 1:length(AV_endpoint_data)
AV_endpoint_train{i} = AV_endpoint_data{i}(1:2:end);%odd trials
AV_endpoint_test{i} = AV_endpoint_data{i}(2:2:end);%even trials
end

%% fit CI model
%sandboxing stuff, for temporary use
%restricted_trials = AV_tar_pairs(:,1) ~= 6 & AV_tar_pairs(:,1) ~= -6;
%
options = optimset('MaxFunEvals',10000,'MaxIter',20000,'PlotFcns',@optimplotfval);
free_params_vector = cell2mat(struct2cell(free_params)); %need to cast params as vector for fminsearch

% run fminsearch to fit parameters
CI_minsearch = @(free_params_vector)get_nll_CI(AV_tar_pairs,AV_endpoint_train,fixed_params,free_params_vector);
[fit_params_vector_CI,min_nll_CI,~,~] = fminsearch(CI_minsearch,free_params_vector,options);

%get nll on test dataset, using fit params
test_nlls.CI = get_nll_CI(AV_tar_pairs,AV_endpoint_test,fixed_params,fit_params_vector_CI);

%put params back in structure format
fit_params_CI = free_params;
names = fieldnames(free_params);
for i=1:length(names)
fit_params_CI.(names{i})= fit_params_vector_CI(i);
end

%% Alternative models
%alternative models for comparing to full CI model
%% fully segregated model, 4 params
free_params_vector_seg = free_params_vector(1:4);
seg_minsearch = @(free_params_vector_seg)get_nll_seg(AV_tar_pairs,AV_endpoint_train,fixed_params,free_params_vector_seg);
[fit_params_vector_seg,min_nll_seg,~,~] = fminsearch(seg_minsearch,free_params_vector_seg,options);

%get nll on test dataset, using fit params
test_nlls.seg = get_nll_seg(AV_tar_pairs,AV_endpoint_test,fixed_params,fit_params_vector_seg);


%% fully integrated model, 4 params
free_params_vector_int = free_params_vector(1:4);
int_minsearch = @(free_params_vector_int)get_nll_int(AV_tar_pairs,AV_endpoint_train,fixed_params,free_params_vector_int);
[fit_params_vector_int,min_nll_int,~,~] = fminsearch(int_minsearch,free_params_vector_int,options);

%get nll on test dataset, using fit params
test_nlls.int = get_nll_int(AV_tar_pairs,AV_endpoint_test,fixed_params,fit_params_vector_int);


%% probability matching instead of posterior reweighted



%% CI model with target locs and sigmas determined by unimodal fits, 2 params
% free_params_vector_fAV = free_params_vector(4:end);
% 
% CI_fAV_minsearch = @(free_params_vector_fAV)get_nll_CI_fixedAV(AV_tar_pairs,AV_endpoint_train,fixed_params,free_params_vector_fAV);
% [fit_params_vector_fAV,min_nll_CI_fAV,~,~] = fminsearch(CI_fAV_minsearch,free_params_vector_fAV,options);
 
%notes: model does not converge, ends up with extremely large prior sig and
%extremely low p_common indicating that these variables are adding very
%little. expectation is that the lack of accuracy is really hurting this
%model and the parameters are trying to account for that somehow. 

%% model comparison

n_obs = length(vertcat(AV_endpoint_test{:}));
n_params.CI = 5;
n_params.seg = 4;
n_params.int = 4;

% put resultsi n table for easier viewing
model_comp_table = get_model_comp_table(test_nlls,n_params,n_obs)

%% generate figures
%plotting demos
if make_plots
plot_dir = 'results\fit_dists';
for i = 1:size(AV_tar_pairs,1)
A_tar = AV_tar_pairs(i,1);
V_tar =  AV_tar_pairs(i,2);
plot_fit_dists(data,A_tar,V_tar,fit_params_CI,fixed_params,-40:40);
saveas(gcf,sprintf('%s\\H01%dA%dV',plot_dir,A_tar,V_tar),'png')
end
end
