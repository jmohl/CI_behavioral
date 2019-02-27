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

global MAXRNG
MAXRNG = 50;
binsize = 1; %size of bins used for responses, in degrees.

% set models to be run
%model list
% [2 1 1] = bayesian CI (2) unity judgement(1) with numerical integration (1)
% [3 1 2] = probabilistic fusion (3) unity (1) analytic solution (2)
% [2 2 1] = bayesian CI (2) localization(2) with numerical integration (1)
% [3 2 2] = probabilistic fusion (3) localization(2) analytic solution (2)
% [2 3 1] = bayesian CI (2) joint fit(3) numerical integration(1)
% [3 3 2] = probabilistic fusion (2) joint fit(3) analytic solution(2)

model_list = {[2,1,1];[3,1,2];[2,2,1];[3,2,2];[2,3,1];[3,3,2]};

%setting fitting procedure options
fitoptions.correct_bias = 1;
fitoptions.load_saved_fits = 0; %load saved fits, if they exist
fitoptions.make_plots = 1;
fitoptions.UBND = [15 15 40 .9 .9]; %upper bounds on theta for grid search
fitoptions.LBND = [1 1 1 .1 0]; %lower bounds on theta
fitoptions.grid_fineness = 3; %number of points per parameter in grid search, remember n points in grid = grid_fineness^n_params;; 
fitoptions.fmin_options = optimset('MaxFunEvals',1000,'MaxIter',500);
fitoptions.eval_range = linspace(-MAXRNG,MAXRNG,MAXRNG*2/binsize + 1); %note can adjust fineness of binning here if wanted. This makes 1 degree bins
fitoptions.eval_midpoints = linspace(-MAXRNG+binsize/2,MAXRNG-binsize/2,length(fitoptions.eval_range)-1);

fitoptions.n_pooled_days = 10; % for using monkey datasets with several days of data

seed = 'default';
%todo:fix this
run_days_separately =0;


if ~exist('results\modelfits', 'dir')
    mkdir('results\modelfits')
end

% load example day, for testing
% raw_data = load('H08_AVD2_2018_08_10_tidy.mat');% load('Juno_AVD2_2017_07_31_2_tidy.mat')
% raw_data = raw_data.tidy_data;
% subject = 'Juno';
%% run model on all subjects
for i= 1:length(subject_list)
    subject = subject_list{i};
    % load data
    if strcmp(subject,'Juno') | strcmp(subject,'Yoko')
        raw_data = load_pool_data(fitoptions.n_pooled_days,subject,seed); %data is pooled across N randomly selected days, yielding a single tidy data table
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
            run_subject
        end
    else
        run_subject
    end
end




