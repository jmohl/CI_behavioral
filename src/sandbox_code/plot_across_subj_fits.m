%% plots comparing subjects model fits
%
% -------------------
% Jeff Mohl
% 9/13/18
% -------------------
%
% Description: compare the model fits across subjects, specifically how
% much better the CI model does for each subject based on AIC and BIC
% compared to integrate and segregate models, as well as the across subject
% average. Plots a horizontal box plot as well as the individual
% performances for each subject (humans) or day (juno)
%
% note: a little hacky with how the box plots are labeled, should probably
% fix at some point but expect this plot to see some changes
%
function plot_across_subj_fits
addpath('results\model_comp')
%load data
Juno_files=dir(fullfile('results\model_comp','*Juno*AVD2*2018*.mat'));
human_files = dir(fullfile('results\model_comp','*H*.mat'));

juno_data = [];
for i = 1:length(Juno_files);
    this_data = load(Juno_files(i).name);
    Juno_IDs{i} = strtok(Juno_files(i).name,'.');
    juno_data = [juno_data;this_data.model_comp_table];
end
human_data = [];
for i = 1:length(human_files);
    this_data = load(human_files(i).name);
    human_IDs{i} = strtok(human_files(i).name,'.');
    human_data = [human_data;this_data.model_comp_table];
end

CI_data_J = juno_data(strcmp(juno_data.Model, 'CI'),:);
CI_data_J.subj = Juno_IDs';
seg_data_J = juno_data(strcmp(juno_data.Model, 'seg'),:);
seg_data_J.subj = Juno_IDs';
int_data_J = juno_data(strcmp(juno_data.Model, 'int'),:);
int_data_J.subj = Juno_IDs';

CI_data_H = human_data(strcmp(human_data.Model, 'CI'),:);
CI_data_H.subj = human_IDs';
seg_data_H = human_data(strcmp(human_data.Model, 'seg'),:);
seg_data_H.subj = human_IDs';
int_data_H = human_data(strcmp(human_data.Model, 'int'),:);
int_data_H.subj = human_IDs';


%make plot for juno and humans
figure()
sc_plot = scatter(seg_data_J.BIC_diff - CI_data_J.BIC_diff,ones(length(Juno_IDs),1),100);
sc_plot.CData = 1:length(Juno_IDs);
sc_plot.MarkerFaceColor = 'flat';
sc_plot.MarkerFaceAlpha = 0.5;
hold on;
sc_plot = scatter(seg_data_H.BIC_diff - CI_data_H.BIC_diff,ones(length(human_IDs),1)*2,100);
sc_plot.CData = 1:length(human_IDs);
sc_plot.MarkerFaceColor = 'flat';
sc_plot.MarkerFaceAlpha = 0.5;
labeled_data = [seg_data_J.BIC_diff - CI_data_J.BIC_diff, ones(height(CI_data_J),1)];
labeled_data = [labeled_data;seg_data_H.BIC_diff - CI_data_H.BIC_diff, ones(height(CI_data_H),1)*2];
boxplot(labeled_data(:,1),labeled_data(:,2), 'orientation', 'horizontal','plotstyle','compact','Colors','k','Labels',{'Juno','Humans'});
xlabel('BIC difference from CI model')
plot([0 0], [0,3],'LineWidth',2,'Color',[0 0 0 .5])
ylim([0 3])
title('Evidence for segregating model vs CI model')
saveas(gcf,'results\model_comp\segVsCI','svg')

%make plot for juno and humans
figure()
sc_plot = scatter(int_data_J.BIC_diff - CI_data_J.BIC_diff,ones(length(Juno_IDs),1),100);
sc_plot.CData = 1:length(Juno_IDs);
sc_plot.MarkerFaceColor = 'flat';
sc_plot.MarkerFaceAlpha = 0.5;
hold on;
sc_plot = scatter(int_data_H.BIC_diff - CI_data_H.BIC_diff,ones(length(human_IDs),1)*2,100);
sc_plot.CData = 1:length(human_IDs);
sc_plot.MarkerFaceColor = 'flat';
sc_plot.MarkerFaceAlpha = 0.5;
labeled_data = [int_data_J.BIC_diff - CI_data_J.BIC_diff, ones(height(CI_data_J),1)];
labeled_data = [labeled_data;int_data_H.BIC_diff - CI_data_H.BIC_diff, ones(height(CI_data_H),1)*2];
boxplot(labeled_data(:,1),labeled_data(:,2), 'orientation', 'horizontal','plotstyle','compact','Colors','k','Labels',{'Juno','Humans'});
plot([0 0], [0,3],'LineWidth',2,'Color',[0 0 0 .5])
ylim([0 3])
title('Evidence for segregating model vs CI model')
xlabel('BIC difference from CI model')
xlim([-200 1500])
saveas(gcf,'results\model_comp\intVsCI','svg')


end
