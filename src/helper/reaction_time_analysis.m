%% Reaction time analysis
%
% -------------------
% Jeff Mohl
% 5/28/20
% -------------------
%
% Description: checking to see if the reaction time differs as a function
% of target separation. Also checking to see if the threshold is smaller
% for single saccades in slower vs faster saccades


%% copied from figure_script. Remove when making into function.
 local_directory = 'D:\GitHub\CI_behavioral\';
 cd(local_directory)

addpath('src', 'src\plotting','results','data','src\lautils');
figpath = 'results\reaction_time_plot';
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

%concatenate human data
H_list = dir('data\H*');
data_H = [];
for fi = 1:length(H_list)
    this_H = load(H_list(fi).name);
    this_H = this_H.tidy_data;
    data_H = vertcat(data_H,this_H);
end
%% plot all reaction times for juno

% data = data_j;
% data=data(data.valid_tr ==1,:);
% 
% %calculate reaction time
% clear reaction_time
% for i = 1:height(data)
%     this_sacs = data(i,:).sac_intervals{:};
%     reaction_time(i) = this_sacs(find(this_sacs(:,1) > data(i,:).go_time-100,1),1)-data(i,:).go_time;
% end
% 
% data.rt = reaction_time';
% data=data(data.rt > 0,:);
% 
% figure; histogram(data.rt) %two peaks indicate that monkey is often trying to anticipate go cue

%% take the slowest 50% of saccades and the fastests 50%, and compute p_single by distance for each
%for each data day, perform this analysis, then compute SEM across days.
%Repeat for each dataset

perc = 50; %what percent of saccades to include in fast and slow categories
subjects = {'MJ','MY','H'};
data = {data_j,data_y,data_H};
for isubj = 1:length(subjects)
subject = subjects{isubj};
this_data = data{isubj};
this_data.sep = this_data.A_tar - this_data.V_tar;
this_data = this_data(abs(this_data.sep) <= 18,:); %don't include the most well separated, because that condition is rare and when subsetting the data here causes problems
this_data = this_data(strcmp(this_data.trial_type,'AV'),:); %only including AV trials for this analysis, because comparing single vs dual saccades
this_data = this_data(this_data.valid_tr ==1,:);

%get reaction times
clear reaction_time
for itr = 1:height(this_data)
    this_sacs = this_data(itr,:).sac_intervals{:};
    reaction_time(itr) = this_sacs(find(this_sacs(:,1) > this_data(itr,:).go_time-100,1),1)-this_data(itr,:).go_time;
end
this_data.rt = reaction_time';
this_data=this_data(this_data.rt > 0,:); %cut anticipatory trials cause those are ambiguous.

%get mean response patterns
g = findgroups(this_data.file_ID);
results_mat_slow = zeros(max(g),length(unique(this_data.sep)));
results_mat_fast = zeros(max(g),length(unique(this_data.sep)));

for ig = 1:max(g)
    
AVdata = this_data(g==ig,:); %select only one day worth of data

slow_data = AVdata(AVdata.rt > prctile(AVdata.rt,100-perc),:);
fast_data = AVdata(AVdata.rt < prctile(AVdata.rt,perc),:);

get_p_common = @(x) sum(x == 1)/ length(x);

[gs,glab] = findgroups(slow_data.sep);
sep_means_slow = splitapply(get_p_common,slow_data.n_sacs,gs);
results_mat_slow(ig,:) = sep_means_slow;
[gs,glab] = findgroups(fast_data.sep);
sep_means_fast = splitapply(get_p_common,fast_data.n_sacs,gs);
results_mat_fast(ig,:) = sep_means_fast;
end

figure
hold on
sem_slow = std(results_mat_slow,1)/sqrt(size(results_mat_slow,1));
sem_fast =  std(results_mat_fast,1)/sqrt(size(results_mat_fast,1));
errorbar(glab, mean(results_mat_slow,1),sem_slow,'r','LineWidth',2)
errorbar(glab, mean(results_mat_fast,1),sem_fast,'g','LineWidth',2) 
title(sprintf('%s, split by rt',subject))
ylabel('p single saccade')
xlabel('target separation')
legend('slow','fast')

saveas(gcf,sprintf('%s\\%s_rt_all',figpath,subject),'svg')
%include 2 way anova for testing purposes
A = reshape(results_mat_fast,[1,size(results_mat_fast,1)*size(results_mat_fast,2)]);
B = reshape(results_mat_slow,[1,size(results_mat_fast,1)*size(results_mat_fast,2)]);
anova_mat = [A',B'];
anova2(anova_mat,size(results_mat_fast,1),'on')
end

