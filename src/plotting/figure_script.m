%% A collection of figures for manuscript
%
% -------------------
% Jeff Mohl
% 7/15/19
% -------------------
%
% Description: generates plots for the CI behavioral manuscript
%

 local_directory = 'C:\Users\jtm47\Documents\Projects\CI_behavioral\';
 cd(local_directory)

addpath('src', 'src\plotting','results','data','src\lautils');
figpath = 'doc\doc_figures';
try
    mkdir(figpath)
end

set(0, 'DefaultLineLineWidth', 1.5);
set(0, 'DefaultLineMarkerSize', 7);
set(0,'DefaultFigurePosition',[25,50,800,800])

%selecting some nice colors for plotting, consistent across all
global model_color aud_color vis_color
model_color = [35/256,155/256,86/256]; %[35/256,155/256,86/256] green
aud_color = [192/256, 57/256, 43/256];
vis_color = [52/256, 152/256, 219/256];

savefiles = 1;

%% load model structures - subsequent plots all rely on fit model results
subject_m = {'Yoko' 'Juno'};
subject_h = {'H02' 'H03' 'H04' 'H05' 'H06' 'H07' 'H08'};

clear models_m models_h models_my models_mj ind;

modelfits_loc = 'modelfits';
%load monkey data
ind = 1;
for subject = subject_m
    m=load(sprintf('results\\%s\\%s_m.mat',modelfits_loc,subject{:}));
    models_m{ind}=m.m;
    ind = ind+1;
    %for each monkey subject, get the individual runs names as well
end
%split up yoko days
file_name = sprintf('/Yoko*AVD2*.mat');
ind_days = dir([sprintf('results\\%s',modelfits_loc) file_name]);
subject_my = {ind_days.name};
ind = 1;
for subject = subject_my
    m=load(sprintf('results\\%s\\%s',modelfits_loc,subject{:}));
    models_my{ind}=m.m;
    ind = ind+1;
end

%split up juno days
file_name = sprintf('/Juno*AVD2*.mat');
ind_days = dir([sprintf('results\\%s',modelfits_loc) file_name]);
subject_mj = {ind_days.name};
ind = 1;
for subject = subject_mj
    m=load(sprintf('results\\%s\\%s',modelfits_loc,subject{:}));
    models_mj{ind}=m.m;
    ind = ind+1;
end

%combined human subjects
ind = 1;
for subject = subject_h
    file_name = sprintf('/%s*AVD2*.mat',subject{:});
    fileid = dir([sprintf('results\\%s',modelfits_loc) file_name]);
    m=load(sprintf('results\\%s\\%s',modelfits_loc,fileid.name));
    models_h{ind}=m.m;
    ind = ind+1;
end

%% Plot Schematics - this runs fairly slow because it is doing numerical integration
plot_schematics(figpath)

%% plot raw data figures - uses raw tidy_data
plot_raw_behavior(figpath,savefiles);

%% get raw data statistics for accuracy
data_j = load('data\Juno_combined.mat');
data_j = data_j.tidy_data;

data_y = load('data\Yoko_combined.mat');
data_y = data_y.tidy_data;

files_H = dir('data\*H*tidy*');
data_H = {};
for ind = 1:length(files_H)
    fname = files_H(ind).name;
    this_data = load(sprintf('data\\%s',fname));
    data_H = vertcat(data_H,this_data.tidy_data);
end

subj = {data_j,data_y,data_H};
subj_name = {'Juno','Yoko','Human'};

for subj_ind = 1:length(subj)
    raw_behavior_stats(subj_ind) = get_raw_behavior_stats(subj{subj_ind},subj_name{subj_ind});
end

for subj_ind = 1:length(subj)
raw_behavior_stats_table(subj_ind,:) = cell2table({raw_behavior_stats(subj_ind).subject,...
    mean(raw_behavior_stats(subj_ind).A_error),...
    mean(raw_behavior_stats(subj_ind).V_error),...
    mean(raw_behavior_stats(subj_ind).std_A),...
    mean(raw_behavior_stats(subj_ind).std_V),...
    mean(raw_behavior_stats(subj_ind).AV_A_error),...
    mean(raw_behavior_stats(subj_ind).AV_V_error),...
    mean(raw_behavior_stats(subj_ind).std_AV_A),...
    mean(raw_behavior_stats(subj_ind).std_AV_V)});
end
raw_behavior_stats_table.Properties.VariableNames = {'Subject','A_error','V_error','std_A','std_V','AV_A_error','AV_V_error','std_AV_A','std_AV_V'};

%% unity judgement plots

% Set demo models to use
% model(1) = CI type: none(0), Bayes (1) or probabilistic fusion (2)
% model(2) = stimulus fusion: none(0), Bayes reweight(1), model selection(2), probabilistic fusion (3)
% model(3) = task/fit type: unity judgement (1), localization (2), joint fit (3), unisensory localization (4)
% model(4) = prior type: naive normal (1), discrete empirical (2), normal mixture empirical (3)
model = [1 1 3 1];

% plot individually for monkeys
% for ind = 1:length(models_m)
%     m=models_m{ind};
%     subject = m.subject;
%     model_ind = ismember(vertcat(m.models{:}),model,'rows');
%     conditions = m.conditions{model_ind};
%     if model(3) == 3 %joint fit models have cell array rather than vectors for these
%         responses = m.responses{model_ind}{1};
%         fit_dist = m.fit_dist{model_ind}{1};
%     else
%         responses = m.responses{model_ind};
%         fit_dist = m.fit_dist{model_ind};
%     end
%     figure
%     plot_unity(conditions,responses,fit_dist);
%     title(sprintf('percent trials reported unity by target separation \n %s model %d%d%d%d',subject, model))
%     if savefiles
%         saveas(gcf,sprintf('%s\\%s_unity_combined',figpath,subject),'svg');
%     end
% end

% plot combined for humans
figure
plot_unity_combined(models_h,model)
if savefiles
    saveas(gcf,sprintf('%s\\human_unity',figpath),'svg');
    saveas(gcf,sprintf('%s\\human_unity',figpath),'png');
end

%plot combined for yoko
figure
plot_unity_combined(models_my,model)
title('% of trials reported unity by target separation - Yoko')
if savefiles
    saveas(gcf,sprintf('%s\\yoko_unity_splitday',figpath),'svg');
    saveas(gcf,sprintf('%s\\yoko_unity_splitday',figpath),'png');
end

% plot combined for juno
figure
plot_unity_combined(models_mj,model)
title('% of trials reported unity by target separation - Juno')
if savefiles
    saveas(gcf,sprintf('%s\\juno_unity_splitday',figpath),'svg');
    saveas(gcf,sprintf('%s\\juno_unity_splitday',figpath),'png');
end


%% localization plots with model fits
% Set demo models to use
% model(1) = CI type: none(0), Bayes (1) or probabilistic fusion (2)
% model(2) = stimulus fusion: none(0), Bayes reweight(1), model selection(2), probabilistic fusion (3)
% model(3) = task/fit type: unity judgement (1), localization (2), joint fit (3), unisensory localization (4)
% model(4) = prior type: naive normal (1), discrete empirical (2), normal mixture empirical (3)
model = [1 1 3 1];
example_conds_h = [13:15];
example_conds_m = [12,14,15];
plot_pred = 1;
%plot for monkey days with pooled model fit

% for ind = 1:length(models_m)
%     m=models_m{ind};
%     subject = m.subject;
%     %generate plot for single subject
%      for ind2 = 1:2:length(example_conds) %for plotting all combos
%     plot_localization(m,model,example_conds([ind2, ind2+1]),plot_pred);
%     set(gcf,'Position',[25,50,1400,500])
%     end
%     if savefiles
%         saveas(gcf,sprintf('%s\\%s_loc_combined',figpath,subject),'svg');
%         saveas(gcf,sprintf('%s\\%s_loc_combined',figpath,subject),'png');
%     end
% end
%if given an array of model fits, will average them together and plot that

% for ind2 = 1:2:length(example_conds_h)  %for plotting all combos
plot_localization(models_h,model,example_conds_h,plot_pred);
set(gcf,'Position',[25,50,1400,300])
% end
if savefiles
    saveas(gcf,sprintf('%s\\humans_loc_combined',figpath),'svg');
    saveas(gcf,sprintf('%s\\humans_loc_combined',figpath),'png');
end

% for ind2 = 1:2:length(example_conds_m)  %for plotting all combos
plot_localization(models_mj,model,example_conds_m,plot_pred);
set(gcf,'Position',[25,50,1400,300])
% end
if savefiles
    saveas(gcf,sprintf('%s\\juno_loc_splitday',figpath),'svg');
    saveas(gcf,sprintf('%s\\juno_loc_splitday',figpath),'png');
end

% for ind2 = 1:2:length(example_conds_m)  %for plotting all combos
plot_localization(models_my,model,example_conds_m,plot_pred);
set(gcf,'Position',[25,50,1400,300])
% end
if savefiles
    saveas(gcf,sprintf('%s\\yoko_loc_splitday',figpath),'svg');
    saveas(gcf,sprintf('%s\\yoko_loc_splitday',figpath),'png');
end

%% Table: comparison of fit parameters
model = [1 1 3 1];
subj_models = {models_h;models_mj;models_my;};
subj_labels = {'Humans','Monkey J', 'Monkey Y'};
theta_labels = {'A_sig','V_sig','prior_sig','p_common','lambda'};
theta_means = zeros(3,5);
theta_sem = zeros(3,5);

for ind = 1:length(subj_models)
    m=subj_models{ind,:};
    subj_thetas = zeros(length(m),5);
    for subj_ind = 1:length(m)
        this_m = m{subj_ind};
        model_ind = ismember(vertcat(this_m.models{:}),model,'rows');
        subj_thetas(subj_ind,:) = this_m.thetas{model_ind};
    end
    theta_means(ind,:) = mean(subj_thetas);
    theta_sem(ind,:) = std(subj_thetas)./sqrt(size(subj_thetas,1));
end

mean_table=array2table(theta_means,'VariableNames',theta_labels,'RowNames',subj_labels)
sem_table=array2table(theta_sem,'VariableNames',theta_labels,'RowNames',subj_labels)

%% ANOVA of percent single saccade by target sep

%% Model comparison using VBA toolbox to get BOR and exceedence probabilities
joint_models = {[1 1 3 1];[1 2 3 1];[1 3 3 1]}; %Bayes, model selection, probabilistic fusion (null)
unity_models = {[1 0 1 1];[2 0 1 1]};
loc_models = {[1 1 2 1];[1 2 2 1];[1 3 2 1]};

models = {unity_models,loc_models,joint_models};
n_params = {[5,5],[5, 5, 5],[5,5,5]};
titles = {'Unity','Localization','Joint'};

figure
fit_results_h = {};
fit_results_m = {};
for ind = 1:length(models)

    model = models{ind};
    n_param = n_params{ind};
    % code for get_model_comp_table returns the nll as well as the AIC and BIC
    
    [AICh,BICh,nll_table_h] = get_model_comp_table(models_h,model,n_param);
    [AICm,BICm,nll_table_m] = get_model_comp_table(models_m,model,n_param);
    %but this is in the wrong direction (subjects x models) so need to
    %transpose
    %here I need the log model evidence, which is equivalent to the BIC or
    %AIC * -0.5
    Lh = -0.5*BICh';
    Lm = -0.5*BICm';
    %Use VBA toolbox to get posterior frequency and protected exceedance
    %input is a KxN array of loglikelihood with K models by N subjects
    options.DisplayWin = 0;
    [~,out_h] = VBA_groupBMC(Lh,options); %this is working, I want to use the posterior.r distribution to report. Also get exceedance probabilities as out.ep and protexted as pxp
    [~,out_m] = VBA_groupBMC(Lm,options);
    fit_results_h{ind} = out_h;
    fit_results_m{ind} = out_m;
    % going to plot protected exceedance probability out.pxp, and model
    % frequency as an error bar plot out.Ef Vf
    subplot(1,length(models),ind)
    hold on
    this_bar = bar([out_h.pxp',out_m.pxp']);
    bar_x_coords = 1:length(out_h.pxp);
    xticks(bar_x_coords);
    bar_x_coords = [bar_x_coords - .14,bar_x_coords + .14];
    errorbar(bar_x_coords, [out_h.Ef',out_m.Ef'], [out_m.Vf',out_m.Vf'], 'k.','MarkerSize',10) %these x values are hard coded, to align the bars and errors, but thats not great.
    plot([0.5 length(model)+.5],[1/length(model), 1/length(model)],'r--')
    %formatting
    if ind == 1
        legend('Humans','Monkeys', 'Posterior Frequency','Chance')
    end
    clear model_names;
    for name_ind = 1:length(model)
        model_names{name_ind} = get_model_names(model{name_ind});
    end
    set(gca,'TickLabelInterpreter','none')
    xticklabels(horzcat(model_names{:}));
    xlabel('Model')
    ylabel('Probability')
    title(titles{ind});
    text(length(model)-.5,.5,sprintf('BOR H: %1.2f\nBOR M: %1.2f',out_h.bor, out_m.bor))
    this_bar(1).FaceColor = [.65 .65 .65];
    this_bar(2).FaceColor = [.85 .85 .85];
end
set(gcf,'Position',[25,50,1400,300])
if savefiles
    saveas(gcf,sprintf('%s\\model_comp',figpath),'png');
    saveas(gcf,sprintf('%s\\model_comp',figpath),'svg');
end


