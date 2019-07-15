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
set(0, 'DefaultLineMarkerSize', 20);
set(0,'DefaultFigurePosition',[25,50,800,800])


%% Plot Schematics - this runs fairly slow because it is doing numerical integration
plot_schematics(figpath)

%% Set demo models to use
% [1 1 3 1] [1 2 3 1] [1 3 3 1]
% model(1) = CI type: none(0), Bayes (1) or probabilistic fusion (2)
% model(2) = stimulus fusion: none(0), Bayes reweight(1), model selection(2), probabilistic fusion (3)
% model(3) = task/fit type: unity judgement (1), localization (2), joint fit (3), unisensory localization (4)
% model(4) = prior type: naive normal (1), discrete empirical (2), normal mixture empirical (3)

subject_m = {'Yoko' 'Juno'};
subject_h = {'H02' 'H03' 'H04' 'H05' 'H06' 'H07' 'H08'};
clear models_m models_h ind;
ind = 1;
for subject = subject_m
    m=load(sprintf('results\\modelfits\\%s_m.mat',subject{:}));
    models_m{ind}=m.m;
    ind = ind+1;
end
ind = 1;
for subject = subject_h
    m=load(sprintf('results\\modelfits\\%s_m.mat',subject{:}));
    models_h{ind}=m.m;
    ind = ind+1;
end

%% Plot behavioral data with model fits

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



%% [DEFUNCT] f1: schematics under different CI models

%going to use the 2d plots for different actual conditions, projected down
%onto the different axes for visualization. Noticed that these
%distributions don't look great at the coarseness I have evaluated them at,
%so going to create a demo pdf with closer spacing and save that for
%plotting

%demo models to use
% [1 1 3 1] [1 2 3 1] [1 3 3 1]
% model(1) = CI type: Bayes (1) or probabilistic fusion (2)
% model(2) = stimulus fusion: Bayes reweight(1), model selection(2), probabilistic fusion (3)
% model(3) = task/fit type: unity judgement (1), localization (2), joint fit (3), unisensory localization (4)
% model(4) = prior type: naive normal (1), discrete empirical (2), normal mixture empirical (3)
m=load(sprintf('results\\modelfits\\%s_m.mat',subject));
m=m.m;

model = [1 1 1 1];
demo_conds = [3,4];

locations = m.fitoptions.eval_midpoints;
model_ind = ismember(vertcat(m.models{:}),model,'rows');

fit_dist = m.fit_dist{model_ind};
if iscell(fit_dist)
    fit_dist = fit_dist{2};
end

for ic = demo_conds
    figure;
    this_dist = squeeze(fit_dist(ic,:,:));
    imagesc(locations,locations,this_dist)
    set(gca, 'YDir', 'normal')
    xlim([-40 40])
    ylim([-40 40])
    xlabel('Auditory Location')
    ylabel('Visual Location')
    title(sprintf('%d Aud, %d Vis pair', m.conditions{model_ind}(ic,:)))
    %     saveas(gcf,sprintf('%s\\%s%dA%dV_pdf',figpath,subject, m.conditions{model_ind}(ic,:)),'png');
    %project into auditory dimension
    I_mat = logical(eye([length(locations),length(locations)]));
    aud_dist = sum(this_dist,1);
    vis_dist = sum(this_dist,2);
    AV_dist = this_dist(I_mat);
    figure
    hold on
    plot(locations,aud_dist,'r')
    plot(locations,vis_dist,'b')
    plot(locations,AV_dist,'k')
    
end
