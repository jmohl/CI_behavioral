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
%todo update this list


%% set initial state
local_directory = 'C:\Users\jtm47\Documents\Projects\CI_behavioral\';
cd(local_directory)
addpath('src', 'src\plotting','results','data');
figpath = 'results\figures';
try
    mkdir(figpath)
end

%% Model and subject to plot
subject = 'Yoko';
model = [1 1 3 1];

m=load(sprintf('results\\modelfits\\%s_m.mat',subject));
m=m.m;
model_ind = ismember(vertcat(m.models{:}),model,'rows');
conditions = m.conditions{model_ind};

%% plot jp1: unity judgement plots grouped by target sep, and across targets
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

    saveas(gcf,sprintf('%s\\%s_unity_combined',figpath,subject),'png');
%% plot jp2: location as distance between A and V saccades for each condition grouped by target sep in the same way.
if model(3) == 3 %joint fit models have cell array rather than vectors for these
    responses = m.responses{model_ind}{2};
    fit_dist = m.fit_dist{model_ind}{2};
else
    responses = m.responses{model_ind};
    fit_dist = m.fit_dist{model_ind};
end

figure
plot_sac_by_sep(conditions,responses,fit_dist);
    saveas(gcf,sprintf('%s\\%s_delt_loc_combined',figpath,subject),'png');
%% plot 1: example 2d joint distribution plots
%load fit data

subject = 'Juno';
m=load(sprintf('results\\modelfits\\%s_m.mat',subject));
m=m.m;
% set desired model and range
model = [1 1 2 1];%keep an eye on this in the future, might change
locations = m.fitoptions.eval_range;
model_ind = ismember(vertcat(m.models{:}),model,'rows');

fit_dist = m.fit_dist{model_ind};
if iscell(fit_dist)
    fit_dist = fit_dist{2};
end
ex_conds = [1:20];
for ic = ex_conds
    figure;
    imagesc(locations,locations,squeeze(fit_dist(ic,:,:)))
    set(gca, 'YDir', 'normal')
    xlim([-40 40])
    ylim([-40 40])
    xlabel('Auditory Location')
    ylabel('Visual Location')
    title(sprintf('%d Aud, %d Vis pair', m.conditions{model_ind}(ic,:)))
%     saveas(gcf,sprintf('%s\\%s%dA%dV_pdf',figpath,subject, m.conditions{model_ind}(ic,:)),'png');
end

%actual data
resp_dist = m.responses{model_ind};
if iscell(resp_dist)
    resp_dist = resp_dist{2};
end
ex_conds = [1:20];
for ic = ex_conds
    figure;
    imagesc(locations,locations,squeeze(resp_dist(ic,:,:)))
    set(gca, 'YDir', 'normal')
    xlim([-40 40])
    ylim([-40 40])
    xlabel('Auditory Location')
    ylabel('Visual Location')
    title(sprintf('%d Aud, %d Vis pair', m.conditions{model_ind}(ic,:)))
%     saveas(gcf,sprintf('%s\\%s%dA%dV_resp',figpath,subject, m.conditions{model_ind}(ic,:)),'png');
end
%% plot 2: unity judgement plots, average across subjects/days

%build unity judgement array of response/fit vectors
subjects_H = {'H02_m.mat' 'H03_m.mat' 'H04_m.mat' 'H05_m.mat' 'H06_m.mat' 'H07_m.mat' 'H08_m.mat'};
J_files = struct2cell(dir('results\modelfits\Juno*'));
subjects_J = J_files(1,:);
Y_files = struct2cell(dir('results\modelfits\Yoko*'));
subjects_Y = Y_files(1,:);
subjects = {subjects_H,subjects_J,subjects_Y};
data_labels = {'Human','J','Y'};
model = [1 0 1];%joint fit model

for data_ind = 1:3
    this_subjects = subjects{data_ind};
    si=1;
    all_resp = zeros(20,length(this_subjects));
    all_fits = all_resp;
    for subject = this_subjects
        m = load(sprintf('results\\modelfits\\%s',subject{:}));
        m=m.m;
        model_ind = ismember(vertcat(m.models{:}),model,'rows');
        if model(3) == 3
            responses = m.responses{model_ind}{1};
            fit_dist = m.fit_dist{model_ind}{1};
        else
            responses = m.responses{model_ind};
            fit_dist = m.fit_dist{model_ind};
        end
        conditions = m.conditions{model_ind};
        all_resp(:,si) = responses(:,1)./sum(responses,2); %just take the percent single (
        %plot_psingle(responses,conditions,fit_dist);
        all_fits(:,si) = fit_dist(:,1);
        si = si +1;
    end
    
    %group values by target separation %TODO
    %sep = conditions(:,1) - conditions(:,2);
    %[sep_g, sep_val] = findgroups(sep);
    
    mean_fit =  mean(all_fits,2); %mean across humans
    fit_CI_int = 1.96 * std(all_fits,0,2)/size(this_subjects,2); % 95% confidence interval
    mean_resp = mean(all_resp,2);
    std_resp = std(all_resp,0,2);
    sem_resp = std_resp / sqrt(length(this_subjects));
    
    % split data by condition, 24 or 6 aud tar, for easier visualization
    A_tars = unique(abs(conditions(:,1)));
    for A_ind = A_tars'
        figure
        hold on;
        pos_inds = conditions(:,1) == A_ind;
        neg_inds = conditions(:,1) == -A_ind;
        v_tar = conditions(pos_inds,2); %using visual target location ,will almost certainly change in future but this provides maximum information
        for subj = 1:size(this_subjects,2)
            %plotting individual responses
            plot(conditions(pos_inds,2),all_resp(pos_inds,subj),'LineWidth',.5,'Color',[.9,.9,.9],'LineStyle','-');
            plot(conditions(neg_inds,2),all_resp(neg_inds,subj),'LineWidth',.5,'Color',[.9,.9,.9],'LineStyle','-');
        end
        p = errorbar(v_tar,mean_resp(pos_inds),sem_resp(pos_inds),'LineWidth',2);
        plot(v_tar,mean_fit(pos_inds),'LineWidth',2,'Color',p.Color,'LineStyle','--');
        %plot(v_tar,mean_fit(pos_inds)+fit_CI_int(pos_inds),'LineWidth',1,'Color',p.Color,'LineStyle','-');
        %plot(v_tar,mean_fit(pos_inds)-fit_CI_int(pos_inds),'LineWidth',1,'Color',p.Color,'LineStyle','-');
        v_tar = conditions(neg_inds,2);
        p = errorbar(v_tar,mean_resp(neg_inds),sem_resp(neg_inds),'LineWidth',2);
        plot(v_tar,mean_fit(neg_inds),'LineWidth',2,'Color',p.Color,'LineStyle','--');
        title(sprintf('Single saccade ratio by tar sep (%s):\\pm %d Aud loc',data_labels{data_ind},A_ind),'Interpreter','tex')
        ylabel('percent single saccade')
        xlabel('Vis target location')
        set(gca,'box','off')
        %set(gcf,'Position',[100,60,1049,895])
        saveas(gcf,sprintf('%s\\psing_%s_%dA_%s',figpath,data_labels{data_ind},A_ind,string(get_model_names(model))),'png');
    end
    
end
%% plot 3: localization plots + models, rescaled and sized for nice figures, separate individuals
subjects = {'Yoko'};%,'Yoko','H05','H03','H08'};
% set desired model and range
model = [1 1 3 1];

for si = 1:length(subjects)
    subject = subjects{si};
    m=load(sprintf('results\\modelfits\\%s_m.mat',subject));
    m=m.m;
    if contains(subject,'H')
        tar_pairs = {[-24,-24];[-24, -18];[-24,-12];[-12,-24];[-12, -18];[-12,-12]}; %tar pairs for humans
    else
        tar_pairs = {[-24,-24];[-24, -12];[-24,-6];[-6,-24];[-6, -18];[-6,-6]};
    end
    
    xlocs = m.fitoptions.eval_midpoints;
    model_ind = ismember(vertcat(m.models{:}),model,'rows');
    if model(3) == 3
        saccades_all = m.responses{model_ind}{2};
        predicted_all = m.fit_dist{model_ind}{2};
    else
        saccades_all = m.responses{model_ind};
        predicted_all = m.fit_dist{model_ind};
    end
    conditions = m.conditions{model_ind};
    
    figure;
    set(gcf,'Position',[0,0,1600,1000])
    for pind = 1:length(tar_pairs)
        subplot(2,3,pind)
        tar_pair = tar_pairs{pind};
        this_ind = ismember(conditions,tar_pair,'rows');
        saccades = saccades_all(this_ind,:,:);
        predicted = predicted_all(this_ind,:,:);
        norm_saccades=saccades/sum(saccades(:));
        norm_predicted = predicted*abs(xlocs(1)-xlocs(2))^2;
        %single saccades will be along the diagonal
        I_mat = logical(eye([length(xlocs),length(xlocs)]));
        sing_sacs = norm_saccades(:,I_mat);
        norm_saccades(:,I_mat) = 0; %remove single saccades from norm sac mat.
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
        legend('Single Sac','A sac','V sac','Model','Location','Best');
        xlabel('endpoint location (degrees)');
        ylabel('% saccades in bin')
        set(gca,'box','off')
        %set bounds that focus on the data instead of the whole range
        %min_tar = min(tar_pair);
        %max_tar = max(tar_pair);
        %xlim([min_tar - 15, max_tar + 15])
        xlim([-30,10])
    end
    saveas(gcf,sprintf('%s\\locex_%s_combined_%s',figpath,subject,string(get_model_names(model))),'png');
end
%% plot 3, take 2, plotting every pair from one side.
subject = 'Juno';
A_tars = [-24 -6];
V_tars = [-24 -18 -12 -6];
A_tars = A_tars * -1; V_tars = V_tars * -1;
% set desired model and range
model = [1 1 3 3];
m=load(sprintf('results\\modelfits\\%s_m.mat',subject));
m=m.m;

xlocs = m.fitoptions.eval_midpoints;
model_ind = ismember(vertcat(m.models{:}),model,'rows');
if model(3) == 3
    saccades_all = m.responses{model_ind}{2};
    predicted_all = m.fit_dist{model_ind}{2};
else
    saccades_all = m.responses{model_ind};
    predicted_all = m.fit_dist{model_ind};
end
conditions = m.conditions{model_ind};

figure;
set(gcf,'Position',[100,60,1600,500])
pind = 1;
for V_tar = V_tars
    for A_ind = A_tars
        subplot(4,2,pind)
        tar_pair = [A_ind,V_tar];
        this_ind = ismember(conditions,tar_pair,'rows');
        saccades = saccades_all(this_ind,:,:);
        predicted = predicted_all(this_ind,:,:);
        norm_saccades=saccades/sum(saccades(:));
        norm_predicted = predicted*abs(xlocs(1)-xlocs(2))^2;
        %single saccades will be along the diagonal
        I_mat = logical(eye([length(xlocs),length(xlocs)]));
        sing_sacs = norm_saccades(:,I_mat);
        norm_saccades(:,I_mat) = 0; %remove single saccades from norm sac mat.
        A_sacs = sum(norm_saccades,2)/2; %sum over all V saccades, divide by 2 to make probabilities sum to 1
        V_sacs = sum(norm_saccades,3)/2;
        
        sac_bar = bar(xlocs,[sing_sacs(:),A_sacs(:),V_sacs(:)],'stacked');%normalizing to probability
        sac_bar(1).FaceColor = [.2 .2 .2];
        sac_bar(2).FaceColor = [1 .5 .5];
        sac_bar(3).FaceColor = [.5 .5 1];
        hold on
        projected_pred = (squeeze(sum(norm_predicted,2)) + squeeze(sum(norm_predicted,3))')/2; %divide by 2 to normalize after marginalizing over both dimesions.
        plot(xlocs,projected_pred,'LineWidth',1.5,'Color',[.5 .5 .5]);
        
        %title(sprintf('%d A, %d V, %s', tar_pair,subject));
        %legend('Single Sac','A sac','V sac','Model','Location','Best');
        xlabel(A_ind);
        ylabel(V_tar)
        set(gca,'box','off')
        %set bounds that focus on the data instead of the whole range
        min_tar = min(tar_pair);
        max_tar = max(tar_pair);
        xlim([min_tar - 15, max_tar + 15])
        %xlim([-30,15])
        pind = pind + 1;
    end
end
saveas(gcf,sprintf('%s\\locex_%s_left_%s',figpath,subject,string(get_model_names(model))),'png');
%% plot 4, splitting up to compare fits in each dimension
subject = 'Yoko';
A_tars = [-24 -6];
V_tars = [-24 -18 -12 -6];
A_tars = A_tars * -1; V_tars = V_tars * -1;
% set desired model and range
model = [1 1 3 3];
m=load(sprintf('results\\modelfits\\%s_m.mat',subject));
m=m.m;

xlocs = m.fitoptions.eval_midpoints;
model_ind = ismember(vertcat(m.models{:}),model,'rows');
if model(3) == 3
    saccades_all = m.responses{model_ind}{2};
    predicted_all = m.fit_dist{model_ind}{2};
else
    saccades_all = m.responses{model_ind};
    predicted_all = m.fit_dist{model_ind};
end
conditions = m.conditions{model_ind};

try
    mkdir(sprintf('%s\\split_by_type',figpath))
end
for V_tar = V_tars
    for A_ind = A_tars
        tar_pair = [A_ind,V_tar];
        this_ind = ismember(conditions,tar_pair,'rows');
        saccades = saccades_all(this_ind,:,:);
        predicted = predicted_all(this_ind,:,:);
        norm_saccades=saccades/sum(saccades(:));
        norm_predicted = predicted*abs(xlocs(1)-xlocs(2))^2; %normalize by bin size
        % why does this some to less than 1? Something seems wrong.
        %single saccades will be along the diagonal
        I_mat = logical(eye([length(xlocs),length(xlocs)]));
        sing_sacs = norm_saccades(:,I_mat);
        norm_saccades(:,I_mat) = 0; %remove single saccades from norm sac mat.
        A_sacs = sum(norm_saccades,2)/2; %sum over all V saccades, divide by 2 to make probabilities sum to 1
        V_sacs = sum(norm_saccades,3)/2;
        
        figure;
        set(gcf,'Position',[100,60,1600,500])
        
        subplot(1,3,1)
        sac_bar = bar(xlocs,sing_sacs(:));
        sac_bar.FaceColor = [.2 .2 .2];
        hold on;
        projected_pred = norm_predicted(:,I_mat);
        plot(xlocs,projected_pred,'LineWidth',1.5,'Color',[.5 .5 .5]);
        title(sprintf(' %d A, %d V',A_ind, V_tar));
        xlim([-15 35])
        
        norm_predicted(:,I_mat) = 0; %remove points from diagonal
        
        subplot(1,3,2)
        sac_bar = bar(xlocs,A_sacs(:));
        sac_bar(1).FaceColor = [1 .5 .5];
        hold on;
        projected_pred = squeeze(sum(norm_predicted,2)/2);%divide by 2 to normalize
        plot(xlocs,projected_pred,'LineWidth',1.5,'Color',[.5 .5 .5]);
        title(sprintf('%d A',A_ind));
        xlim([-15 35])
              
        subplot(1,3,3)
        sac_bar = bar(xlocs,V_sacs(:));
        sac_bar(1).FaceColor = [.5 .5 1];
        hold on
        projected_pred = squeeze(sum(norm_predicted,3)/2);
        plot(xlocs,projected_pred,'LineWidth',1.5,'Color',[.5 .5 .5]);
        title(sprintf('%d V',V_tar));
        xlim([-15 35])
        
        saveas(gcf,sprintf('%s\\split_by_type\\%s_%d%d',figpath,subject,A_ind,V_tar),'png');
    end
end
%% plot 5 perform model comparison and make plots

%n params is 5 for all models except the probabilistic fusion model in the
%unity judgement case
subject_list = {'Juno' 'Yoko' 'H02' 'H03' 'H04' 'H05' 'H06' 'H07' 'H08'};
n_models = 3;
aic = zeros(length(subject_list),n_models);% number of models needed
bic = aic;
models =[1 1 3; 1 2 3; 1 3 3];
nparams = [5, 5, 5];
for i= 1:length(subject_list)
    m=load(sprintf('results\\modelfits\\%s_m.mat',subject_list{i}));
    m=m.m;
    model_ind = ismember(vertcat(m.models{:}),models,'rows');
    nll=m.nll(model_ind);
    %sum nll across k folds for each model, if relevant
    if m.fitoptions.cross_validate
        nll=[cell2mat(nll{1,1}),cell2mat(nll{1,2}),cell2mat(nll{1,3})];
        nll = sum(nll,1);
        nobs = sum(m.responses{model_ind(1)}{1,1}(:)); %maybe a bug here TODO
    else
        nll = cell2mat(nll);
        nobs = sum(sum(m.responses{1},[]));
    end
    [aic(i,:),bic(i,:)] = aicbic(-nll,nparams,nobs);
end

%get model names using model values
model_names = get_model_names(models);

aic_table = array2table(aic);
aic_table.Properties.VariableNames = strcat(model_names,'_AIC');
bic_table =   array2table(bic);
bic_table.Properties.VariableNames = strcat(model_names,'_BIC');
mc_table = [aic_table,bic_table];
mc_table.Properties.RowNames = subject_list;

%JM todo fix this
dif_bic(:,1) = bic(:,1) - bic(:,2); %bayesian - model selection model, BIC dif, joint fit
dif_bic(:,2) = bic(:,1) - bic(:,3); %bayesian - probabilistic fusion model
dif_bic = array2table(dif_bic);
dif_bic.Properties.VariableNames = {'MS','PF'};
dif_bic.Properties.RowNames = subject_list;

% plot results
figure;
sc_plot = scatter(dif_bic{:,1},ones(length(subject_list),1),100);
sc_plot.CData = 1:length(subject_list);
sc_plot.MarkerFaceColor = 'flat';
sc_plot.MarkerFaceAlpha = 0.5;
hold on
sc_plot = scatter(dif_bic{:,2},ones(length(subject_list),1)*2,100);
sc_plot.CData = 1:length(subject_list);
sc_plot.MarkerFaceColor = 'flat';
sc_plot.MarkerFaceAlpha = 0.5;
labeled_data = [dif_bic{:,1}, ones(length(subject_list),1)];
labeled_data = [labeled_data;dif_bic{:,2},ones(length(subject_list),1)*2];
boxplot(labeled_data(:,1),labeled_data(:,2), 'orientation', 'horizontal','plotstyle','compact','Colors','k','Labels',{'MS','PF'});
xlabel('BIC difference (BB - model)')
plot([0 0], [0,3],'LineWidth',2,'Color',[0 0 0 .5])
ylim([0 3])
title('Evidence for Bayes/Bayes model vs alternatives')
saveas(gcf,sprintf('%s\\bic_dif',figpath),'png');
%% table 1 - table of model fit parameters

subject_list = {'Juno' 'Yoko'};% 'H02' 'H03' 'H04' 'H05' 'H06' 'H07' 'H08'};
model =[1 1 3];
n_params= 5;
params = zeros(length(subject_list),n_params);
params_sd = params;

for i= 1:length(subject_list)
    m=load(sprintf('results\\modelfits\\%s_m.mat',subject_list{i}));%This code is meant to be run on the cross validated data
    m=m.m;
    model_ind = ismember(vertcat(m.models{:}),model,'rows');
    thetas =m.thetas(model_ind);
    if m.fitoptions.cross_validate    %sum nll across k folds for each model
        thetas=cellfun(@cell2mat,thetas,'UniformOutput',0);
        mean_thetas = cellfun(@mean,thetas,'UniformOutput',0);
        std_thetas = cellfun(@std,thetas,'UniformOutput',0);
        params(i,:) = mean_thetas{:};
        params_sd(i,:) = std_thetas{:}; %not actually sure this is what I want. Might want the standard deviation of the mean values rather than the std of the estimates across folds.
    else
        params(i,:) = thetas{:};
        params_sd(i,:) = std(thetas{:});
    end
end

param_table = array2table(params);
param_table.Properties.VariableNames ={'V_sig','A_sig','p_sig','p_common','lambda'};
param_table.Properties.RowNames = subject_list;

param_sd_table = array2table(params_sd);
param_sd_table.Properties.VariableNames ={'V_sig','A_sig','p_sig','p_common','lambda'};
param_sd_table.Properties.RowNames = subject_list;

human_table = param_table(3:end,:);
grpstats(human_table,[],{'mean','std'})
%% unisensory localization

subject = 'Juno';
m=load(sprintf('results\\modelfits\\%s_m.mat',subject));
m=m.m;
model = [0 0 4 3];

xlocs = m.fitoptions.eval_midpoints;
model_ind = ismember(vertcat(m.models{:}),model,'rows');
saccades_V = m.responses{model_ind}{2};
predicted_V = m.fit_dist{model_ind}{2};
saccades_A = m.responses{model_ind}{1};
predicted_A = m.fit_dist{model_ind}{1};
conditions_A = m.conditions{model_ind}{1};
conditions_V = m.conditions{model_ind}{2};

figure
for A_ind = 1:length(conditions_A)
    subplot(2,2,A_ind)
    saccades = saccades_A(A_ind,:);
    predicted = predicted_A(A_ind,:);
    norm_saccades=saccades/sum(saccades(:));
    
    plot(xlocs,norm_saccades,'r','LineWidth',1.5)
    hold on
    plot(xlocs,predicted,'LineWidth',1.5,'Color',[.5 .5 .5]);
    plot([conditions_A(A_ind), conditions_A(A_ind)],[0,.05],'k--','LineWidth',1.5);
    xlabel('Location (deg)')
    ylabel('Percent of saccades')
    title(sprintf('%d Auditory',conditions_A(A_ind)))
    set(gca,'box','off')
end
set(gcf,'Position',[0,0,1400,1400])
saveas(gcf,sprintf('%s\\%s_unimodal_A_%d%d%d%d',figpath,subject,model),'png');

figure
for V_ind = 1:length(conditions_V)
    subplot(5,2,V_ind)
    saccades = saccades_V(V_ind,:);
    predicted = predicted_V(V_ind,:);
    norm_saccades=saccades/sum(saccades(:));
    
    plot(xlocs,norm_saccades,'b','LineWidth',1.5)
    hold on
    plot(xlocs,predicted,'LineWidth',1.5,'Color',[.5 .5 .5]);
    plot([conditions_V(V_ind), conditions_V(V_ind)],[0,.2],'k--','LineWidth',1.5);
    xlabel('Location (deg)')
    ylabel('Percent of saccades')
    title(sprintf('%d Visual',conditions_V(V_ind)))
    set(gca,'box','off')
end
set(gcf,'Position',[0,0,1400,1400])

saveas(gcf,sprintf('%s\\%s_unimodal_V_%d%d%d%d',figpath,subject,model),'png');

