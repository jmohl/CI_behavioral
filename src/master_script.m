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
local_directory = 'C:\Users\jtm47\Documents\Projects\CI_behavioral\';
cd(local_directory)
addpath('data','src','src\lautils', 'src\plotting');

subject_list = {'Juno' 'Yoko' 'H02' 'H03' 'H04' 'H05' 'H06' 'H07' 'H08'};


% set models to be run

%% Select models to run

%new model list
% [1 0 1] = bayesian CI (1) no combination (0) unity judgement task(1)
% [2 0 1] = probabilistic fusion (2) no combination (0) unity (1)
% [1 1 2] = bayesian CI (1) with reweighting (1) localization task(2)
% [2 1 2] = probabilistic fusion (2) with reweighting (1) localization(2) 
% [1 1 3] = bayesian CI (1) reweighting(1) joint fit(3) 
% [1 2 3] = bayesian CI (1) model selection(2) joint fit(3)
% [1 3 3] = bayesian CI (1) probabilistic fusion(3) joint fit(3)
% [0 0 4] = no CI (0), no fusion (0), unimodal localization (4);
model_list = {[1 0 1]; [2 0 1]; [1 1 2]; [2 1 2]; [1 1 3]; [1 2 3]; [1 3 3]; [0 0 4]};
%setting fitting procedure options
global fitoptions MAXRNG
MAXRNG = 50;
fitoptions.binsize = 1; %size of bins used for responses, in degrees.
fitoptions.load_saved_fits = 1; %load saved fits, if they exist
fitoptions.make_plots = 0;
fitoptions.n_iterations = 5; %set option to repeat fminsearch for n times
fitoptions.UBND = [6 15 40 .9 .25]; %upper bounds on theta for grid search
fitoptions.LBND = [.5 1 5 .1 .01]; %lower bounds on theta for grid search
fitoptions.grid_fineness = 4; %number of points per parameter in grid search, remember n points in grid = grid_fineness^n_params;; 
fitoptions.fmin_options = optimset('MaxFunEvals',1000,'MaxIter',1000, 'TolFun', 1e-2, 'TolX',1e-2);
fitoptions.eval_range = linspace(-MAXRNG,MAXRNG,MAXRNG*2/fitoptions.binsize + 1);
fitoptions.eval_midpoints = linspace(-MAXRNG+fitoptions.binsize/2,MAXRNG-fitoptions.binsize/2,length(fitoptions.eval_range)-1);

fitoptions.cross_validate = 0;
fitoptions.kfolds = 5;

fitoptions.dynamic_bins = 0; %experimental
%todo:fix this
run_days_separately = 0;

if ~exist('results\modelfits', 'dir')
    mkdir('results\modelfits')
end

% load example day, for testing
% subject_list = {'Juno', 'Yoko'};
% model_list = {[0 0 4]}; %testing unisensory localization under same framework.

%% run model on all subjects
for i= 1:length(subject_list)
    subject = subject_list{i};
    % load data
    if strcmp(subject,'Juno') | strcmp(subject,'Yoko')
        this_file = dir(sprintf('data\\*%s_combined_30*',subject));
        raw_data=load(this_file.name);
        raw_data = raw_data.tidy_data;
    else
        this_file = dir(sprintf('data\\*%s*',subject));
        load(this_file.name);
        raw_data = tidy_data;
    end
    %TODO: fix this deadline inspired hack
    if run_days_separately
        comb_data = raw_data;
        days_list = unique(raw_data.file_ID);
        for this_day=1:length(days_list)
            raw_data = comb_data(strcmp(comb_data.file_ID, days_list{this_day}),:);
            subject = days_list{this_day};
            run_subject(subject,raw_data,model_list)
        end
    else
        run_subject(subject,raw_data,model_list)
    end
end




