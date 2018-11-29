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
close all;
local_directory = 'C:\Users\jtm47\Documents\Projects\CI_behavioral\';
CI_opts.n_pooled_days = 5; % for using monkey datasets with several days of data
seed = 'default';
CI_opts.make_plots = 1;
CI_opts.correct_bias = 1;
CI_opts.k_folds = 10; 
subject_list = {'Juno' 'Yoko' 'H02' 'H03' 'H04' 'H05' 'H06' 'H07' 'H08'};

%todo:fix this
run_days_separately =0;

%adding paths
cd(local_directory)
addpath('src','results','data')

% load initial parameters
load('ini_params.mat')
%fminsearch options
fmin_options = optimset('MaxFunEvals',20000,'MaxIter',40000);

%% run model on all subjects
for i= 1:length(subject_list)
    subject = subject_list{i};
    % load data
    if strcmp(subject,'Juno') | strcmp(subject,'Yoko')
        raw_data = load_pool_data(CI_opts.n_pooled_days,subject,seed); %data is pooled across N randomly selected days, yielding a single tidy data table
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
    %plot distributions of saccades as well as predicted distributions under CI model fit 
    if CI_opts.make_plots
        plot_diagnostics(fixed_params,data,subject)
        set(0,'DefaultFigureVisible','off');
        plot_dir = sprintf('results\\fit_dists\\%s',subject);
        mkdir(plot_dir);
        for j = 1:size(fixed_params.AV_pairs,1)
            A_tar = fixed_params.AV_pairs(j,1);
            V_tar =  fixed_params.AV_pairs(j,2);
            plot_fit_dists(data,[A_tar V_tar],'CI',fit_params('CI',:),fixed_params,-40:40);
            saveas(gcf,sprintf('%s\\%s%dA%dV',plot_dir,subject,A_tar,V_tar),'png')
        end
        set(0,'DefaultFigureVisible','on');
    end
end

%% generate figures
% % plotting demos, this section
% if CI_opts.make_plots
%     % diagnostic plots for this subject, characterizing unimodal bias
%     plot_diagnostics(fixed_params,data,subject)
%     % plot distributions of saccades as well as predicted distributions under CI model fit 
%     set(0,'DefaultFigureVisible','off');
%     plot_dir = sprintf('results\\fit_dists\\%s',subject);
%     mkdir(plot_dir);
%     for i = 1:size(AV_tar_pairs,1)
%         A_tar = AV_tar_pairs(i,1);
%         V_tar =  AV_tar_pairs(i,2);
%         this_A_data = data.A{A_tars == A_tar};
%         this_V_data = data.V{V_tars == V_tar};
%         this_AV_data = data.AV{i};
%         plot_fit_dists(this_A_data,this_V_data,this_AV_data,A_tar,V_tar,fit_params_CI,fixed_params,-40:40);
%         saveas(gcf,sprintf('%s\\%s%dA%dV',plot_dir,subject,A_tar,V_tar),'png')
%     end
%     set(0,'DefaultFigureVisible','on');
% end
