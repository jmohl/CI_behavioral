%% A collection of figures for manuscript
%
% -------------------
% Jeff Mohl
% 6/14/19
% -------------------
%
% Description: generates following plots

local_directory = 'C:\Users\jtm47\Documents\Projects\CI_behavioral\';
cd(local_directory)
addpath('src', 'src\plotting','results','data');
figpath = 'results\figures';
try
    mkdir(figpath)
end

%% f1: schematics under different CI models

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

model = [1 1 3 1];
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
