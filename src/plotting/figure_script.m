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

set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultLineMarkerSize', 7);
set(0,'DefaultFigurePosition',[25,50,800,800])

savefiles = 1;

%% Plot Schematics - this runs fairly slow because it is doing numerical integration
 plot_schematics(figpath)

%% plot raw data figures - uses raw tidy_data
 plot_raw_behavior(figpath,savefiles)

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

%% unity judgement plots
%TODO: 1: use joint model fit 2: update to work with cross validated?

% Set demo models to use
% model(1) = CI type: none(0), Bayes (1) or probabilistic fusion (2)
% model(2) = stimulus fusion: none(0), Bayes reweight(1), model selection(2), probabilistic fusion (3)
% model(3) = task/fit type: unity judgement (1), localization (2), joint fit (3), unisensory localization (4)
% model(4) = prior type: naive normal (1), discrete empirical (2), normal mixture empirical (3)
model = [1 0 1 1];

% plot individually for monkeys
for ind = 1:length(models_m)
    m=models_m{ind};
    subject = m.subject;
    model_ind = ismember(vertcat(m.models{:}),model,'rows');
    conditions = m.conditions{model_ind};
    if model(3) == 3 %joint fit models have cell array rather than vectors for these
        responses = m.responses{model_ind}{1};
        fit_dist = m.fit_dist{model_ind}{1};
    else
        responses = m.responses{model_ind};
        fit_dist = m.fit_dist{model_ind};
    end
    figure
    plot_unity(conditions,responses,fit_dist);
    title(sprintf('percent trials reported unity by target separation \n %s model %d%d%d%d',subject, model))
    if savefiles
        saveas(gcf,sprintf('%s\\%s_unity_combined',figpath,subject),'svg');
    end
end

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
model = [1 1 2 1];
example_conds_h = [11 13];
example_conds_m = [11 12];
plot_pred = 1;
%plot for monkey days with pooled model fit

% for ind = 1:length(models_m)
%     m=models_m{ind};
%     subject = m.subject;
%     %generate plot for single subject
%      for ind2 = 1:2:length(example_conds) %for plotting all combos
%     plot_localization(m,model,example_conds([ind2, ind2+1]),plot_pred);
%     set(gcf,'Position',[25,50,1300,500])
%     end
%     if savefiles
%         saveas(gcf,sprintf('%s\\%s_loc_combined',figpath,subject),'svg');
%         saveas(gcf,sprintf('%s\\%s_loc_combined',figpath,subject),'png');
%     end
% end
%if given an array of model fits, will average them together and plot that

for ind2 = 1:2:length(example_conds_h)  %for plotting all combos
    plot_localization(models_h,model,example_conds_h([ind2, ind2+1]),plot_pred);
    set(gcf,'Position',[25,50,1300,500])
end
if savefiles
    saveas(gcf,sprintf('%s\\humans_loc_combined',figpath),'svg');
    saveas(gcf,sprintf('%s\\humans_loc_combined',figpath),'png');
end

for ind2 = 1:2:length(example_conds_m)  %for plotting all combos
    plot_localization(models_mj,model,example_conds_m([ind2, ind2+1]),plot_pred);
    set(gcf,'Position',[25,50,1300,500])
end
if savefiles
    saveas(gcf,sprintf('%s\\juno_loc_splitday',figpath),'svg');
    saveas(gcf,sprintf('%s\\juno_loc_splitday',figpath),'png');
end

for ind2 = 1:2:length(example_conds_m)  %for plotting all combos
    plot_localization(models_my,model,example_conds_m([ind2, ind2+1]),plot_pred);
    set(gcf,'Position',[25,50,1300,500])
end
if savefiles
    saveas(gcf,sprintf('%s\\yoko_loc_splitday',figpath),'svg');
    saveas(gcf,sprintf('%s\\yoko_loc_splitday',figpath),'png');
end
%% Condensed localization plot

model = [1 1 2 1];
true_loc = 0; %option to use true target locations or relative locations (from unimodal saccades) for specifying bias.
plot_condensed_loc(models_mj,model,true_loc);
title('Juno')
if savefiles
    saveas(gcf,sprintf('%s\\juno_condensed_bias',figpath),'png');
    saveas(gcf,sprintf('%s\\juno_condensed_bias',figpath),'svg');
end

plot_condensed_loc(models_my,model,true_loc);
title('Yoko')
if savefiles
    saveas(gcf,sprintf('%s\\yoko_condensed_bias',figpath),'png');
    saveas(gcf,sprintf('%s\\yoko_condensed_bias',figpath),'svg');
end

plot_condensed_loc(models_h,model,true_loc);
title('Human')
if savefiles
    saveas(gcf,sprintf('%s\\HU_condensed_bias',figpath),'png');
    saveas(gcf,sprintf('%s\\HU_condensed_bias',figpath),'svg');
end


%% unimodal localization plot, not working yet
% example_conds = [2 4];
% plot_localization(models_h,[0 0 4 1],example_conds,plot_pred);

%% AIC BIC table
%which models to compare
models = {[1 1 3 1];[1 2 3 1];[1 3 3 1]}; %Bayes, model selection, probabilistic fusion (null) 
%comparing models on the localization component only, as the 
n_params = [5, 5, 5];

[AIC_mj,BIC_mj] = get_model_comp_table(models_mj,models,n_params);
[AIC_my,BIC_my] = get_model_comp_table(models_my,models,n_params);
[AIC_h,BIC_h] = get_model_comp_table(models_h,models,n_params);

%get relative to non-CI model (probabilistic fusion)
AIC_mj_rel = AIC_mj - AIC_mj(:,3);
AIC_my_rel = AIC_my - AIC_my(:,3);
AIC_h_rel = AIC_h - AIC_h(:,3);

BIC_mj_rel = BIC_mj - BIC_mj(:,3);
BIC_my_rel = BIC_my - BIC_my(:,3);
BIC_h_rel = BIC_h - BIC_h(:,3);

%get means + std for all
%rows are [monkey_j;monkey_y;humans]
%columns are [models 1:3]
AIC_mean = [mean(AIC_mj_rel,1);mean(AIC_my_rel,1);mean(AIC_h_rel,1)];
AIC_sem = [std(AIC_mj_rel,1)/sqrt(length(models_mj));std(AIC_my_rel,1)/sqrt(length(models_my));std(AIC_h_rel,1)/sqrt(length(models_h))];
BIC_mean = [mean(BIC_mj_rel,1);mean(BIC_my_rel,1);mean(BIC_h_rel,1)];
BIC_sem = [std(BIC_mj_rel,1)/sqrt(length(models_mj));std(BIC_my_rel,1)/sqrt(length(models_my));std(BIC_h_rel,1)/sqrt(length(models_h))];

for ind = 1:length(models)
model_names{ind} = get_model_names(models{ind});
end
model_comp_AIC = array2table(AIC_mean],'RowNames',{'mj','my','h'},'VariableNames',horzcat(model_names{:}))
model_comp_BIC = array2table(BIC_mean],'RowNames',{'mj','my','h'},'VariableNames',horzcat(model_names{:}))

    
