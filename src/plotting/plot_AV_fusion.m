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
AV_color = [.5 .2 .5];

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

data_h8 = load('data\H08_AVD2_2018_08_10_tidy.mat'); %4_17, 4_25 good examples but have non-0 fix
data_h8 = data_h8.tidy_data;
%% main plots

this_data = data_y;
this_subj = 'Yoko';
this_type = 'sep';%coi, sep, or all
this_data = this_data(this_data.valid_tr==1,:);
%might also want to adjust so AV sacs are only single sacs?
this_data = this_data(this_data.n_sacs == 1,:);

% pairs = [-6 -6;-6 -12;-24 -24; -24 -18;6 6; 6 12; 24 24; 24 18];
all_pairs = table2array(unique(this_data(strcmp(this_data.trial_type,'AV'),{'A_tar','V_tar'}))); % to get all target pairs
if strfind(this_subj,'H')
    %for humans
coi_pairs = [-24 -24; -12 -12; 12 12; 24 24];
sep_pairs = [-24 -18; -12 -6; 12 6; 24 18];
else
coi_pairs = [-24 -24; -6 -6; 6 6; 24 24];
sep_pairs = [-24 -18; -6 -12; 6 12; 24 18];
end

if strcmp('coi', this_type)
    pairs = coi_pairs;
elseif strcmp('sep', this_type)
    pairs = sep_pairs;
else
    pairs=all_pairs;
end
close all

means_fig = figure(21);
title(sprintf('Subject: %s, tars: %s',this_subj,this_type))
hold on;
set(means_fig,'Position',[500,500,600,200])

bins = [-30:30];
stats_arr = zeros(2,3,length(pairs));

for ind = 1:length(pairs)
this_pair = pairs(ind,:);

this_AV_sac = vertcat(this_data(this_data.A_tar == this_pair(1) & this_data.V_tar == this_pair(2),:).valid_endpoints{:});
this_V_sac = vertcat(this_data(this_data.V_tar == this_pair(2) & strcmp(this_data.trial_type,'V'),:).valid_endpoints{:});
this_A_sac = vertcat(this_data(this_data.A_tar == this_pair(1) & strcmp(this_data.trial_type,'A'),:).valid_endpoints{:});

%plot histograms for vis, comparing AV, A, and V. Only x coord

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

%plot mean and STD, also compare with optimal, for a bunch of target pairs
figure(means_fig)
errorbar(V_mean,1,V_std,'horizontal','LineStyle','none','LineWidth',2,'Color',vis_color);
errorbar(AV_mean,2,AV_std,'horizontal','LineStyle','none','LineWidth',2,'Color',AV_color);
errorbar(A_mean,3,A_std,'horizontal','LineStyle','none','LineWidth',2,'Color',aud_color);
scatter([V_mean,AV_mean,A_mean],[1:3],20,[vis_color; AV_color;aud_color],'filled')

plot([this_pair(1) this_pair(1)],[0,4],'--', 'Color',aud_color)
plot([this_pair(2) this_pair(2)],[0,4],'--','Color', vis_color)

xlim([-30,30])
ylim([0,4])
yticks([1:3])
yticklabels({'Vsac','AVsac','Asac'});
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

%% comparing AV singles with AV doubles
this_data = data_h8;
this_subj = 'H08';
this_type = 'sep';%coi, sep, or all
this_data = this_data(this_data.valid_tr==1,:);
%might also want to adjust so AV sacs are only single sacs?
this_dual_data = this_data(this_data.n_sacs > 1,:);
this_data = this_data(this_data.n_sacs == 1,:);

% pairs = [-6 -6;-6 -12;-24 -24; -24 -18;6 6; 6 12; 24 24; 24 18];
all_pairs = table2array(unique(this_data(strcmp(this_data.trial_type,'AV'),{'A_tar','V_tar'}))); % to get all target pairs
if strfind(this_subj,'H')
    %for humans
coi_pairs = [-24 -24; -12 -12; 12 12; 24 24];
sep_pairs = [-24 -18; -12 -18; 12 18; 24 18];
else
coi_pairs = [-24 -24; -6 -6; 6 6; 24 24];
sep_pairs = [-24 -12; -6 -18; 6 18; 24 12];
end

if strcmp('coi', this_type)
    pairs = coi_pairs;
elseif strcmp('sep', this_type)
    pairs = sep_pairs;
else
    pairs=all_pairs;
end
close all
bins = [-30:30];

for ind = 1:length(pairs)
this_pair = pairs(ind,:);

this_AV_sac = vertcat(this_data(this_data.A_tar == this_pair(1) & this_data.V_tar == this_pair(2),:).valid_endpoints{:});
this_AV_sep_A = vertcat(this_dual_data(this_dual_data.A_tar == this_pair(1) & this_dual_data.V_tar == this_pair(2),:).A_endpoints{:});
this_AV_sep_V = vertcat(this_dual_data(this_dual_data.A_tar == this_pair(1) & this_dual_data.V_tar == this_pair(2),:).V_endpoints{:});

%plot histograms for vis, comparing AV, A, and V. Only x coord

figure
histogram(this_AV_sep_V(:,1),bins,'Normalization','probability')
hold on
histogram(this_AV_sep_A(:,1),bins,'Normalization','probability')
histogram(this_AV_sac(:,1),bins,'Normalization','probability')
legend('V sac', 'A sac', 'AV sac')
title(sprintf('%dA %dV',this_pair))

%stats
AV_mean = mean(this_AV_sac(:,1));
AV_std = std(this_AV_sac(:,1));
V_mean = mean(this_AV_sep_V(:,1));
V_std = std(this_AV_sep_V(:,1));
A_mean = mean(this_AV_sep_A(:,1));
A_std = std(this_AV_sep_A(:,1));
figure
title(sprintf('Subject: %s, tars: %s',this_subj,this_type))
hold on;
set(gcf,'Position',[500,500,600,200])
%plot mean and STD, also compare with optimal, for a bunch of target pairs
errorbar(V_mean,1,V_std,'horizontal','LineStyle','none','LineWidth',2,'Color',vis_color);
errorbar(AV_mean,2,AV_std,'horizontal','LineStyle','none','LineWidth',2,'Color',AV_color);
errorbar(A_mean,3,A_std,'horizontal','LineStyle','none','LineWidth',2,'Color',aud_color);
scatter([V_mean,AV_mean,A_mean],[1:3],20,[vis_color; AV_color;aud_color],'filled')

plot([this_pair(1) this_pair(1)],[0,4],'--', 'Color',aud_color)
plot([this_pair(2) this_pair(2)],[0,4],'--','Color', vis_color)

xlim([-30,30])
ylim([0,4])
yticks([1:3])
yticklabels({'Vsac','AVsac','Asac'});

end

