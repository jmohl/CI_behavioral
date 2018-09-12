%% script for running a single daaset through the analysis
%
% -------------------
% Jeff Mohl
% 9/6/18
% -------------------
%
% Description: this script will allow running of each step of the
% behavioral CI project, from loading data, fitting model, and generating
% figures, for a single subject. It is meant to be run out of
% master_script and so at least the hard coded section of that script is
% assumed to have been run.

%note: H01 has some problems and it is best to exclude for now, since I
%will probably have to rerun that subject in any case

%todo: put this in function form to reduce workspace clutter


%% load data
if strcmp(subject,'Juno')
    raw_data = load_pool_data(CI_opts.n_pooled_days,seed); %data is pooled across N randomly selected days, yielding a single tidy data table
else
    this_file = dir(sprintf('data\\*%s*',subject));
    load(this_file.name);
    raw_data = tidy_data;
end

raw_data.valid_endpoints = get_response_endpoints(raw_data,0,100)';
%omit trials without valid endpoints, almost always occurs when trial was aborted
data = raw_data(~cellfun('isempty',raw_data.valid_endpoints),:);

%pull out specific locations used in this dataset
fixed_params.V_tars = sortrows(unique(data.V_tar(~isnan(data.V_tar))));
fixed_params.A_tars = sortrows(unique(data.A_tar(~isnan(data.A_tar))));
fixed_params.AV_pairs = sortrows(unique([data.A_tar(~isnan(data.V_tar)& ~isnan(data.A_tar)),data.V_tar(~isnan(data.V_tar)&~isnan(data.A_tar))],'rows'));

%% correct eye tracker calibration
if CI_opts.correct_bias
[data] = get_bias_corrected_data(data);
end

%% split data into train and test sets
[AV_train,AV_test] = get_train_test(data,CI_opts.k_folds,fixed_params.AV_pairs);

%% fit CI model
param_vector = cell2mat(struct2cell(free_params)); %need to cast params as vector for fminsearch

% run fminsearch to fit parameters
for this_fold = 1:CI_opts.k_folds
CI_minsearch = @(free_param_vector)get_nll_CI(AV_train(this_fold,:),fixed_params,free_param_vector);
[fit_results.CI_fit(this_fold,:),fit_results.CI_nll,~,~] = fminsearch(CI_minsearch,param_vector,fmin_options);

%get nll on test dataset, using fit params
test_nlls.CI(this_fold) = get_nll_CI(AV_test(this_fold,:),fixed_params,fit_results.CI_fit(this_fold,:));
end

fit_params_CI = free_params;
names = fieldnames(free_params);
for i=1:length(names)
fit_params_CI.(names{i})= mean(fit_results.CI_fit(:,i));
end
%% Alternative models
%alternative models for comparing to full CI model

% fully segregated model, 4 params
param_vector = param_vector(1:4);
for this_fold = 1:CI_opts.k_folds
seg_minsearch = @(param_vector)get_nll_seg(AV_train(this_fold,:),fixed_params,param_vector);
[fit_results.seg_fit,fit_results.seg_nll,~,~] = fminsearch(seg_minsearch,param_vector,fmin_options);

%get nll on test dataset, using fit params
test_nlls.seg(this_fold) = get_nll_seg(AV_test(this_fold,:),fixed_params,fit_results.seg_fit);


% fully integrated model, 4 params
param_vector = param_vector(1:4);
int_minsearch = @(param_vector)get_nll_int(AV_train(this_fold,:),fixed_params,param_vector);
[fit_results.int_fit,fit_results.int_nll,~,~] = fminsearch(int_minsearch,param_vector,fmin_options);

%get nll on test dataset, using fit params
test_nlls.int(this_fold) = get_nll_int(AV_test(this_fold,:),fixed_params,fit_results.int_fit);
end

%% probability matching instead of posterior reweighted
% todo

%% CI model with target locs and sigmas determined by unimodal fits, 2 params
% get single modality estimates (todo)
% [fixed_params.V_mu,fixed_params.V_sig] = get_unimodal_est(data.V);
% [fixed_params.A_mu,fixed_params.A_sig] = get_unimodal_est(data.A);
% free_params_vector_fAV = free_params_vector(4:end);
% 
% CI_fAV_minsearch = @(free_params_vector_fAV)get_nll_CI_fixedAV(AV_tar_pairs,AV_endpoint_train,fixed_params,free_params_vector_fAV);
% [fit_params_vector_fAV,min_nll_CI_fAV,~,~] = fminsearch(CI_fAV_minsearch,free_params_vector_fAV,options);
 
%notes: model does not converge, ends up with extremely large prior sig and
%extremely low p_common indicating that these variables are adding very
%little. expectation is that the lack of accuracy is really hurting this
%model and the parameters are trying to account for that somehow. 

%% model comparison
n_params.CI = 5;
n_params.seg = 4;
n_params.int = 4;

% put results in table for easier viewing
subject
model_comp_table = get_model_comp_table(test_nlls,n_params,length(vertcat(AV_test{:})))
save(sprintf('results\\model_comp\\mc_%s.mat',subject),'model_comp_table');
writetable(model_comp_table,sprintf('results\\model_comp\\mc_%s.csv',subject));

