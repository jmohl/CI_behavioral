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


%% Plot Schematics - this runs fairly slow because it is doing numerical integration
plot_schematics(figpath)

%% load model structures.
subject_m = {'Yoko' 'Juno'};
subject_h = {'H02' 'H03' 'H04' 'H05' 'H06' 'H07' 'H08'};

clear models_m models_h models_my models_mj ind;

modelfits_loc = 'modelfits';
%load monkey data
% ind = 1;
% for subject = subject_m
%     m=load(sprintf('results\\%s\\%s_m.mat',modelfits_loc,subject{:}));
%     models_m{ind}=m.m; 
%     ind = ind+1;
%     %for each monkey subject, get the individual runs names as well 
% 
% end
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
    ind_days = dir([sprintf('results\\%s',modelfits_loc) file_name]);
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
    saveas(gcf,sprintf('%s\\%s_unity_combined',figpath,subject),'svg');
end

% plot combined for humans
figure
plot_unity_combined(models_h,model)
saveas(gcf,sprintf('%s\\human_unity',figpath),'svg');
saveas(gcf,sprintf('%s\\human_unity',figpath),'png');


%plot combined for yoko
figure
plot_unity_combined(models_my,model)
title('% of trials reported unity by target separation - Yoko')
saveas(gcf,sprintf('%s\\yoko_unity_splitday',figpath),'svg');
saveas(gcf,sprintf('%s\\yoko_unity_splitday',figpath),'png');


% plot combined for juno
figure
plot_unity_combined(models_mj,model)
title('% of trials reported unity by target separation - Juno')
saveas(gcf,sprintf('%s\\juno_unity_splitday',figpath),'svg');
saveas(gcf,sprintf('%s\\juno_unity_splitday',figpath),'png');


%% localization plots with model fits
% Set demo models to use
% model(1) = CI type: none(0), Bayes (1) or probabilistic fusion (2)
% model(2) = stimulus fusion: none(0), Bayes reweight(1), model selection(2), probabilistic fusion (3)
% model(3) = task/fit type: unity judgement (1), localization (2), joint fit (3), unisensory localization (4)
% model(4) = prior type: naive normal (1), discrete empirical (2), normal mixture empirical (3)
model = [1 1 2 1];
example_conds = [2 4]; 
plot_pred = 1;

%plot for monkey days with pooled model fit
for ind = 1:length(models_m)
    m=models_m{ind};
    subject = m.subject;
    %generate plot for single subject
    plot_localization(m,model,example_conds,plot_pred);
    set(gcf,'Position',[25,50,1300,500])
%     saveas(gcf,sprintf('%s\\%s_loc_combined',figpath,subject),'svg');
%     saveas(gcf,sprintf('%s\\%s_loc_combined',figpath,subject),'png');
end
%if given an array of model fits, will average them together and plot that
for ind = [5,15]
    example_conds = [5, 15];
plot_localization(models_h,model,example_conds,plot_pred);
set(gcf,'Position',[25,50,1300,500])
% saveas(gcf,sprintf('%s\\humans_loc_combined',figpath),'svg');
% saveas(gcf,sprintf('%s\\humans_loc_combined',figpath),'png');
% 
plot_localization(models_mj,model,example_conds,plot_pred);
set(gcf,'Position',[25,50,1300,500])
% saveas(gcf,sprintf('%s\\juno_loc_splitday',figpath),'svg');
% saveas(gcf,sprintf('%s\\juno_loc_splitday',figpath),'png');

%some issue here where one of the days apparently doesn't have enough data
%for one of the uni conditions. So that's annoying.
plot_localization(models_my,model,example_conds,plot_pred);
set(gcf,'Position',[25,50,1300,500])
% saveas(gcf,sprintf('%s\\yoko_loc_splitday',figpath),'svg');
% saveas(gcf,sprintf('%s\\yoko_loc_splitday',figpath),'png');
end
%% Condensed localization plot

model = [1 1 2 1];
true_loc = 0; %option to use true target locations or relative locations (from unimodal saccades) for specifying bias.
plot_condensed_loc(models_mj,model,true_loc);
title('Juno')
saveas(gcf,sprintf('%s\\juno_condensed_bias',figpath),'png');

plot_condensed_loc(models_my,model,true_loc);
title('Yoko')
saveas(gcf,sprintf('%s\\yoko_condensed_bias',figpath),'png');

plot_condensed_loc(models_h,model,true_loc);
title('Human')
saveas(gcf,sprintf('%s\\HU_condensed_bias',figpath),'png');


%% unimodal localization plot
example_conds = [2 4]; 
plot_localization(models_h,[0 0 4 1],example_conds,plot_pred);

