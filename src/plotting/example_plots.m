%% A collection of figures for presentation, generated from fit data
%
% -------------------
% Jeff Mohl
% 2/18/19
% -------------------
%
% Description: generates following plots
% plot 1: example joint distribution for the two saccade case
% plot 2: unity judgements pooled across subjects (or days)
% plot 3: specific localization fits, rescaled and sized for nice viewing.


%% set initial state
local_directory = 'C:\Users\jtm47\Documents\Projects\CI_behavioral\';
cd(local_directory)
addpath('src', 'src\plotting','results');
figpath = 'results\ex_figs';

%% plot 1: example 2d joint distribution plots
%load fit data
m=load('results\modelfits\Juno_m.mat');
m=m.m;
% set desired model and range
model = [2 2 1];%keep an eye on this in the future, might change
locations = m.fitoptions.eval_range;
model_ind = ismember(vertcat(m.models{:}),model,'rows');

fit_dist = m.fit_dist{model_ind};
ex_conds = [1,5];
for ic = ex_conds
    figure;
    imagesc(locations,locations,squeeze(fit_dist(ic,:,:)))
    set(gca, 'YDir', 'normal')
    xlim([-40 40])
    xlabel('Auditory Location')
    ylabel('Visual Location')
    title(sprintf('%d Aud, %d Vis pair', m.conditions{model_ind}(ic,:)))
    saveas(gcf,sprintf('%s\\%dA%dV',figpath, m.conditions{model_ind}(ic,:)),'png');
end



%% plot 2: unity judgement plots, average across subjects/days

%build unity judgement array of response/fit vectors
subjects_H = {'H02_m.mat' 'H03_m.mat' 'H04_m.mat' 'H05_m.mat' 'H06_m.mat' 'H07_m.mat' 'H08_m.mat'};
J_files = struct2cell(dir('results\modelfits\Juno_AVD2*'));
subjects_J = J_files(1,:);
Y_files = struct2cell(dir('results\modelfits\Yoko_AVD2*'));
subjects_Y = Y_files(1,:);
subjects = {subjects_H,subjects_J,subjects_Y};
data_labels = {'Human','J','Y'};

for data_ind = 1:3
    this_subjects = subjects{data_ind};
    si=1;
    all_resp = zeros(20,length(this_subjects));
    all_fits = all_resp;
for subject = this_subjects
    m = load(sprintf('results\\modelfits\\%s',subject{:}));
    m=m.m;
    model = [2 1 1];%unity judgement model
    model_ind = ismember(vertcat(m.models{:}),model,'rows');
    fit_dist = m.fit_dist{model_ind};
    responses = m.responses{model_ind};
    conditions = m.conditions{model_ind};
    all_resp(:,si) = responses(:,1)./sum(responses,2); %just take the percent single (
    %plot_psingle(responses,conditions,fit_dist);
    all_fits(:,si) = fit_dist(:,1);
    si = si +1;
end

mean_fit =  mean(all_fits,2); % TODO confidence interval instead at some point?
mean_resp = mean(all_resp,2);
std_resp = std(all_resp,0,2);
sem_resp = std_resp / sqrt(length(this_subjects));

% split data by condition, 24 or 6 aud tar, for easier visualization
A_tars = unique(abs(conditions(:,1)));
for A_tar = A_tars'
    figure
    pos_inds = conditions(:,1) == A_tar;
    neg_inds = conditions(:,1) == -A_tar;
    v_tar = conditions(pos_inds,2); %using visual target location ,will almost certainly change in future but this provides maximum information
    p = errorbar(v_tar,mean_resp(pos_inds),sem_resp(pos_inds),'LineWidth',2);
    hold on;
    plot(v_tar,mean_fit(pos_inds),'LineWidth',2,'Color',p.Color,'LineStyle','--');
    
    v_tar = conditions(neg_inds,2);
    p = errorbar(v_tar,mean_resp(neg_inds),sem_resp(neg_inds),'LineWidth',2);
    plot(v_tar,mean_fit(neg_inds),'LineWidth',2,'Color',p.Color,'LineStyle','--');
    title(sprintf('Single saccade ratio by tar sep (%s):\\pm %d Aud loc',data_labels{data_ind},A_tar),'Interpreter','tex')
    ylabel('percent single saccade')
    xlabel('Vis target location')
    set(gca,'box','off')
    %set(gcf,'Position',[100,60,1049,895])
    saveas(gcf,sprintf('%s\\psing_%s_%dA',figpath,data_labels{data_ind},A_tar),'png');
end
end


%% plot 3: localization plots + models, rescaled and sized for nice figures
subject = 'H03';
tar_pairs = {[-12,-12];[-12, -18];[-12,-24]};

m=load(sprintf('results\\modelfits\\%s_m.mat',subject));
m=m.m;
% set desired model and range
model = [2 3 1];
xlocs = m.fitoptions.eval_midpoints;
model_ind = ismember(vertcat(m.models{:}),model,'rows');
if model(2) == 3
saccades_all = m.responses{model_ind}{2};
predicted_all = m.fit_dist{model_ind}{2};
else
saccades_all = m.responses{model_ind};
predicted_all = m.fit_dist{model_ind};
end
conditions = m.conditions{model_ind};

figure;
set(gcf,'Position',[100,60,1600,500])
for pind = 1:length(tar_pairs)
    subplot(1,3,pind)
    tar_pair = tar_pairs{pind};
    this_ind = ismember(conditions,tar_pair,'rows');
    saccades = saccades_all(this_ind,:,:);
    predicted = predicted_all(this_ind,:,:);
    norm_saccades=saccades/sum(saccades(:));
    norm_predicted = predicted*abs(xlocs(1)-xlocs(2))^2;
    %single saccades will be along the diagonal
    I_mat = logical(eye([length(xlocs),length(xlocs)]));
    sing_sacs = norm_saccades(:,I_mat);
    norm_saccades(:,I_mat) = 0; %remove saccades that are AV from counts
    A_sacs = sum(norm_saccades,2)/2; %divide by 2 to make probabilities sum to 1
    V_sacs = sum(norm_saccades,3)/2;
    sac_bar = bar(xlocs,[sing_sacs(:),A_sacs(:),V_sacs(:)],'stacked');%normalizing to probability
    sac_bar(1).FaceColor = [.2 .2 .2];
    sac_bar(2).FaceColor = [1 .5 .5];
    sac_bar(3).FaceColor = [.5 .5 1];
    hold on
    projected_pred = (squeeze(sum(norm_predicted,2)) + squeeze(sum(norm_predicted,3))')/2; %divide by 2 to normalize
    plot(xlocs,projected_pred,'LineWidth',1.5,'Color',[.5 .5 .5]); %scaling by bin width so that the probability is matched correctly
    
    title(sprintf('%d A, %d V, %s', tar_pair,subject));
    legend('Single Sac','A sac','V sac','Modeled','Location','Best');
    xlabel('endpoint location (degrees)');
    ylabel('% saccades in bin')
    set(gca,'box','off')
    %set bounds that focus on the data instead of the whole range
    %min_tar = min(tar_pair);
    %max_tar = max(tar_pair);
    %xlim([min_tar - 15, max_tar + 15])
    xlim([-30,10])
end
    saveas(gcf,sprintf('%s\\loc\\locex_%s_combined%d%d%d',figpath,subject,model),'png');

%% plot 4 perform model comparison and make plots

%n params is 5 for all models except the probabilistic fusion model in the
%unity judgement case
AIC_table = cell2table(subject_list);
BIC_table = cell2table(subject_list);
aic = zeros(length(subject_list),6);% number of models needed
bic = aic;
model_array = [];

for i= 1:length(subject_list)
    m=load(sprintf('results\\modelfits\\%s_m.mat',subject_list{i}));
    m=m.m;
    models(i,:) = cellfun(@mat2str,m.models,'UniformOutput',0);
    nll = horzcat(m.nll{:});
    nparams = repmat(5,1,length(models));
    nparams(strcmp(models, '[3 1 2]')) =1; % this model only has 1 param.
    nobs = sum(m.responses{1}(:));
    [aic(i,:),bic(i,:)] = aicbic(-nll,nparams,nobs);
end
    










