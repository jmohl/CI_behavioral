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
close all; clear all;
local_directory = 'D:\GitHub\CI_behavioral';
cd(local_directory)
addpath('data','src','src\lautils', 'src\plotting');

subject_list = {'Juno' 'Yoko' 'H02' 'H03' 'H04' 'H05' 'H06' 'H07' 'H08'};
run_days_separately = 1;

%% Select models to run

% model values
% model(1) = CI type: Bayes (1) or probabilistic fusion (2)
% model(2) = stimulus fusion: Bayes reweight(1), model selection(2), probabilistic fusion (3)
% model(3) = task/fit type: unity judgement (1), localization (2), joint fit (3), unisensory localization (4)
% model(4) = prior type: naive normal (1), discrete empirical (2), normal mixture empirical (3)
% load example day, for testing
%model_list = {[0 0 4 1];[0 0 4 3];[1 0 1 1];[1 0 1 3];[2 0 1 1]; [1 1 2 1];[1 1 2 3]; [1 2 2 1]; [1 2 2 3]; [1 1 3 1]; [1 1 3 3]; [1 2 3 1];[1 2 3 3]};
%reduced model list, using minimal viable models;
model_list = {[0 0 4 1];[1 0 1 1];[1 1 2 1];[1 1 3 1];[2 0 1 1];[1 2 2 1];[1 2 3 1];[1 3 2 1];[1 3 3 1]};

%setting fitting procedure options
global fitoptions MAXRNG
MAXRNG = 50;
fitoptions.binsize = 1; %size of bins used for responses, in degrees.
fitoptions.load_saved_fits = 1; %load saved fits, if they exist
fitoptions.n_iterations = 1; %set option to repeat fminsearch for n times
fitoptions.parameter_names = {'A_sig','V_sig','prior_sig','p_common','lambda_uni','prior_mu'};
fitoptions.grid_fineness = 3; %number of points per parameter in grid search, remember n points in grid = grid_fineness^n_params;; 
fitoptions.fmin_options = optimset('MaxFunEvals',10000,'MaxIter',10000, 'TolFun', 1e-2, 'TolX',1e-2);
fitoptions.eval_range = linspace(-MAXRNG,MAXRNG,MAXRNG*2/fitoptions.binsize + 1);
fitoptions.eval_midpoints = linspace(-MAXRNG+fitoptions.binsize/2,MAXRNG-fitoptions.binsize/2,length(fitoptions.eval_range)-1);

fitoptions.cross_validate = 0;
fitoptions.kfolds = 5;

fitoptions.strict_filter=0; %experimental, only applies to localization component, implemented in get_prepro_data
fitoptions.use_uni_means= 0; %experimental 7/23/19 all of these have been evaluated and found to not make much of a difference, so will probably exclude for simplicity
fitoptions.dynamic_bins = 0; %experimental

if ~exist('results\modelfits', 'dir')
    mkdir('results\modelfits')
end

%% testing subset
%  model_list = {[1 3 2 1];[0 0 4 1];[1 1 2 1];};
%  subject_list = {'Juno' 'Yoko'};
%  
 %% run model on all subjects
for i= 1:length(subject_list)
    subject = subject_list{i};
    % load data
    if strcmp(subject,'Juno') | strcmp(subject,'Yoko')
        this_file = dir(sprintf('data\\*%s_combined*',subject));
        raw_data=load(this_file.name);
        raw_data = raw_data.tidy_data;
    else
        this_file = dir(sprintf('data\\*%s*',subject));
        load(this_file.name);
        raw_data = tidy_data;
    end
    run_subject(subject,raw_data,model_list)
    % if want to run the monkey days separately, for better direct
    % comparison with humans
    if run_days_separately
        comb_data = raw_data;
        days_list = unique(raw_data.file_ID);
        for this_day=1:length(days_list)
            raw_data = comb_data(strcmp(comb_data.file_ID, days_list{this_day}),:);
            subject = days_list{this_day};
            run_subject(subject,raw_data,model_list)
        end
    end
end




