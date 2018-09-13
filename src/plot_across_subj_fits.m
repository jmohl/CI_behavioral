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
% average
%
function plot_across_subj_fits
addpath('results\model_comp')
%load data
all_files=dir(fullfile('results\model_comp','*.mat'));

all_data = [];
for i = 1:length(all_files);
    this_data = load(all_files(i).name);
    subj_ids{i} = strtok(all_files(i).name,'.');
    all_data = [all_data;this_data.model_comp_table];
end

CI_data = all_data(strcmp(all_data.Model, 'CI'),:);
CI_data.subj = subj_ids';
seg_data = all_data(strcmp(all_data.Model, 'seg'),:);
seg_data.subj = subj_ids';
int_data = all_data(strcmp(all_data.Model, 'int'),:);
int_data.subj = subj_ids';

%make plot
figure()
sc_plot = scatter(ones(length(subj_ids),1),seg_data.BIC_diff - CI_data.BIC_diff);
sc_plot.CData = 1:9;
sc_plot.MarkerFaceColor = 'flat';
sc_plot.MarkerFaceAlpha = 0.5;
hold on;
b_plot = bar(mean(seg_data.BIC_diff - CI_data.BIC_diff));
b_plot.EdgeColor = 'k';
b_plot.LineWidth = 2;
b_plot.FaceColor = 'none';
xlim([0 2]);


end
