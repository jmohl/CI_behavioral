%% AV saccade plots
%
% -------------------
% Jeff Mohl
% 5/22/20
% -------------------
%
% Description: generates plots to evaluate the location of AV single
% saccades vs AV double saccades and vs A and V single saccades. Purpose of
% these plots is to evaluate whether AV single are more consistent with
% fusion (i.e. occuring intermediate of A and V) or if they could plausibly
% be explained by simply ignoring one of the two targets.


%% copied from figure_script. Remove when making into function.
 local_directory = 'D:\GitHub\CI_behavioral\';
 cd(local_directory)

addpath('src', 'src\plotting','results','data','src\lautils');
figpath = 'results\AV_fuse_plot';
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
%% load data
data_j = load('data\Juno_combined.mat');
data_j = data_j.tidy_data;

data_y = load('data\Yoko_combined.mat');
data_y = data_y.tidy_data;

data_j1 = load('data\Juno_AVD2_2018_04_11_tidy.mat'); %3_19_2, 4_11
data_j1 = data_j1.tidy_data;

data_y1 = load('data\Yoko_AVD2_2019_04_10_tidy.mat'); %4_17, 4_25 good examples but have non-0 fix
data_y1 = data_y1.tidy_data;

%%

this_data = data_j;
this_data = this_data(this_data.valid_tr==1,:);
%might also want to adjust so AV sacs are only single sacs?
this_data = this_data(this_data.n_sacs == 1,:);
bins = [-30:30];

pairs = [-6 -6;-6 -12;-24 -24; -24 -18;6 6; 6 12; 24 24; 24 18];
pairs = table2array(unique(this_data(strcmp(this_data.trial_type,'AV'),{'A_tar','V_tar'}))); % to get all target pairs

stats_arr = zeros(2,3,length(pairs));
for ind = 1:length(pairs)
this_pair = pairs(ind,:);


this_AV_sac = vertcat(this_data(this_data.A_tar == this_pair(1) & this_data.V_tar == this_pair(2),:).valid_endpoints{:});
this_V_sac = vertcat(this_data(this_data.V_tar == this_pair(2) & strcmp(this_data.trial_type,'V'),:).valid_endpoints{:});
this_A_sac = vertcat(this_data(this_data.A_tar == this_pair(1) & strcmp(this_data.trial_type,'A'),:).valid_endpoints{:});

%plot histograms for vis, comparing AV, A, and V. IOnly x coord

figure
histogram(this_V_sac(:,1),bins,'Normalization','probability')
hold on
histogram(this_A_sac(:,1),bins,'Normalization','probability')
histogram(this_AV_sac(:,1),bins,'Normalization','probability')
legend('V sac', 'A sac', 'AV sac')
title(sprintf('%dA %dV',this_pair))

%stats
AV_mean = mean(this_AV_sac(:,1));
AV_std = std(this_AV_sac(:,1));
V_mean = mean(this_V_sac(:,1));
V_std = std(this_V_sac(:,1));
A_mean = mean(this_A_sac(:,1));
A_std = std(this_A_sac(:,1));

stats_arr(:,:,ind) = [AV_mean, V_mean, A_mean;AV_std,V_std,A_std];
end

%% looking to see if the difference from V target changes as a function of target sep
%if aud is on the right (positive disp) expect V target shift to right (up
%in plot)
AV_V_diff = squeeze(stats_arr(1,1,:) - stats_arr(1,2,:));
tar_sep = pairs(:,1) - pairs(:,2);

figure
plot(tar_sep,AV_V_diff,'k.')
xlim([-15 15]);

%that does seem to happen, but effect seems weak. Maybe plot I want to show
%for this is doing 1 sided ttest between AV and V distributions.
