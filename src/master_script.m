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
% - cross validation instead of train/test datasets
% - summarize results of cross validation for each subject
% - plot fit parameters for examples on each subject


%% hardcoded parameters
close all;
local_directory = 'C:\Users\jtm47\Documents\Projects\CI_behavioral\';
CI_opts.n_pooled_days = 15; % for using monkey datasets with several days of data
seed = 'default';
CI_opts.make_plots = 0;
CI_opts.correct_bias = 1;
subject_list = {'Juno' 'H02' 'H03' 'H04' 'H05' 'H06' 'H07' 'H08'};

%adding paths
cd(local_directory)
addpath('src','results','data')

% load initial parameters
load('ini_params.mat')
%fminsearch options
fmin_options = optimset('MaxFunEvals',10000,'MaxIter',20000);

%% run model on all subjects
for i= 1:length(subject_list)
    subject = subject_list{i};
    run_subject
end

%% generate figures
% plotting demos, this section
if CI_opts.make_plots
    % diagnostic plots for this subject, characterizing unimodal bias
    plot_diagnostics(fixed_params,subject)
    % plot distributions of saccades as well as predicted distributions under CI model fit 
    set(0,'DefaultFigureVisible','off');
    plot_dir = sprintf('results\\fit_dists\\%s',subject);
    mkdir(plot_dir);
    for i = 1:size(AV_tar_pairs,1)
        A_tar = AV_tar_pairs(i,1);
        V_tar =  AV_tar_pairs(i,2);
        this_A_data = data.A{A_tars == A_tar};
        this_V_data = data.V{V_tars == V_tar};
        this_AV_data = data.AV{i};
        plot_fit_dists(this_A_data,this_V_data,this_AV_data,A_tar,V_tar,fit_params_CI,fixed_params,-40:40);
        saveas(gcf,sprintf('%s\\%s%dA%dV',plot_dir,subject,A_tar,V_tar),'png')
    end
    set(0,'DefaultFigureVisible','on');
end
