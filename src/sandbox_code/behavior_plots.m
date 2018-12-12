%% Demo plots script
%
% -------------------
% Jeff Mohl
% 12/10/18
% -------------------
%
% Description: this script will perform some simple behavioral analytics
% and generate plots for those, as well as plots comparing the fit
% distributions with the actual behavior.

%% hardcoded params
CI_opts.n_pooled_days = 6; % for using monkey datasets with several days of data
seed = 'default';
CI_opts.make_plots = 1;
CI_opts.correct_bias = 1;
CI_opts.k_folds = 10;
human_list = {'H02' 'H03' 'H04' 'H05' 'H06' 'H07' 'H08'};
subject_list = human_list;

addpath('src','results','data','results\fit_params')
run_days_separately = 1;

if run_days_separately
   file_name = sprintf('/*Juno*AVD2*2018*.mat');
   %get the days for which I have saved parameters
   data_dir = dir(['results/fit_params' file_name]);
   for ii = 1:length(data_dir)-1 %removing the weird outlier day, which is the last one
    juno_list{ii} = sprintf('%s_tidy.mat',data_dir(ii).name(strfind(data_dir(ii).name,'Juno'):end-4));
   end
   subject_list = [subject_list,juno_list];
    
   file_name = sprintf('/Yoko*AVD2*2018*.mat');
   data_dir = dir(['data' file_name]);
        rand_days = randsample(length(data_dir),CI_opts.n_pooled_days);
        subject_list = [subject_list,{data_dir(rand_days).name}];
else
    subject_list = [subject_list,{'Yoko','Juno'}];

end

%% generate perc_saccade structure, for plotting percent single saccade by condition
% todo, this code is very slow
if exist('perc_sing_sac') clear 'perc_sing_sac'; end

for ii= 1:length(subject_list)
    subject = subject_list{ii};
    % load data
    if strcmp(subject,'Juno') | strcmp(subject,'Yoko')
        raw_data = load_pool_data(CI_opts.n_pooled_days,subject,seed); %data is pooled across N randomly selected days, yielding a single tidy data table
    else
        this_file = dir(sprintf('data\\*%s*',subject));
        tidy_data = load(this_file.name);
        if isfield(tidy_data,'this_tidy')
            tidy_data = tidy_data.this_tidy; % because of bug in some of the individually made data files
        else
            tidy_data = tidy_data.tidy_data;
        end
        raw_data = tidy_data;
    end
    
    % data cleaning
    raw_data.valid_endpoints = get_response_endpoints(raw_data,0,100)';
    %omit trials without valid endpoints, almost always occurs when trial was aborted
    data = raw_data(~cellfun('isempty',raw_data.valid_endpoints),:);
    
    % correct eye tracker calibration
    if CI_opts.correct_bias
        [data] = get_bias_corrected_data(data);
    end
    if ~exist('perc_sing_sac')
        perc_sing_sac = get_perc_sing_sac(data);
    else
        perc_sing_sac = [perc_sing_sac;get_perc_sing_sac(data)];
    end
            
end

%% make plots
% make plotting directory

try
    mkdir('results\behavior')
end

%get mean and std by target separation, across subjects (humans first)
human_pss = perc_sing_sac(~cellfun(@isempty,strfind(perc_sing_sac.ID,'H')),:);
juno_pss = perc_sing_sac(~cellfun(@isempty,strfind(perc_sing_sac.ID,'Juno')),:);
yoko_pss = perc_sing_sac(~cellfun(@isempty,strfind(perc_sing_sac.ID,'Yoko')),:);

%make plots for humans
[hg,nvp] = findgroups(human_pss(:,'disp'));

h_means = splitapply(@mean,human_pss.psingle,hg);
h_std = splitapply(@std, human_pss.psingle,hg);

h_len = splitapply(@length,human_pss.psingle,hg);
h_sem = h_std./sqrt(h_len);

figure(1)
errorbar(nvp.disp,h_means,h_sem,'LineWidth',2)
title('Single saccade ratio by target separation - human')
ylabel('percent single saccade')
xlabel('degrees of target separation (vis from aud)')
saveas(gcf,'results\behavior\human_psing','png')

%make plots for Juno
[jg,nvp] = findgroups(juno_pss(:,'disp'));

h_means = splitapply(@mean,juno_pss.psingle,jg);
h_std = splitapply(@std, juno_pss.psingle,jg);
h_len = splitapply(@length,juno_pss.psingle,jg);
h_sem = h_std./sqrt(h_len);

figure(2)
errorbar(nvp.disp,h_means,h_sem,'LineWidth',2)
title('Single saccade ratio by target separation - Juno')
ylabel('percent single saccade')
xlabel('degrees of target separation (vis from aud)')
saveas(gcf,'results\behavior\Juno_psing','png')

% make plots again, this time splitting by aud location
[jg,nvp] = findgroups(juno_pss(:,{'disp' 'A_tar'}));

h_means = splitapply(@mean,juno_pss.psingle,jg);
h_std = splitapply(@std, juno_pss.psingle,jg);
h_len = splitapply(@length,juno_pss.psingle,jg);
h_sem = h_std./sqrt(h_len);
figure(22)
hold on
A_tars = unique(nvp.A_tar);
for ii = 1:length(A_tars)
    this_tar = nvp.A_tar == A_tars(ii);
errorbar(nvp(this_tar,:).disp,h_means(this_tar),h_sem(this_tar),'LineWidth',2)

end
legend('-24','-6','6','24')
title('Single saccade ratio by target separation - Juno')
ylabel('percent single saccade')
xlabel('degrees of target separation (vis from aud)')
saveas(gcf,'results\behavior\Juno_psing_split','png')


%make plots for yoko
[yg,nvp] = findgroups(yoko_pss(:,'disp'));

h_means = splitapply(@mean,yoko_pss.psingle,yg);
h_std = splitapply(@std, yoko_pss.psingle,yg);
h_len = splitapply(@length,yoko_pss.psingle,yg);
h_sem = h_std./sqrt(h_len);

figure(3)
errorbar(nvp.disp,h_means,h_sem,'LineWidth',2)
title('Single saccade ratio by target separation - Yoko')
ylabel('percent single saccade')
xlabel('degrees of target separation (vis from aud)')
saveas(gcf,'results\behavior\Yoko_psing','png')

% make plots again, this time splitting by aud location
[yg,nvp] = findgroups(yoko_pss(:,{'disp' 'A_tar'}));

h_means = splitapply(@mean,yoko_pss.psingle,yg);
h_std = splitapply(@std, yoko_pss.psingle,yg);
h_len = splitapply(@length,yoko_pss.psingle,yg);
h_sem = h_std./sqrt(h_len);
figure(33)
hold on
A_tars = unique(nvp.A_tar);
for ii = 1:length(A_tars)
    this_tar = nvp.A_tar == A_tars(ii);
errorbar(nvp(this_tar,:).disp,h_means(this_tar),h_sem(this_tar),'LineWidth',2)

end
legend('-24','-6','6','24')
title('Single saccade ratio by target separation - yoko')
ylabel('percent single saccade')
xlabel('degrees of target separation (vis from aud)')
saveas(gcf,'results\behavior\Yuno_psing_split','png')
%% Get parameter values for each saved subject
%get juno params
for ii = 1:length(juno_list)
    this_params = load(sprintf('fit_params_%s.mat',juno_list{ii}(1:end-9)));
    this_params = this_params.fit_params;
    j_params(ii,:) = this_params('CI',:);
end

%get median param values, because outliers
med_j_params = median(j_params{:,:});

for ii = 1:length(human_list)
    this_params = load(sprintf('fit_params_%s',human_list{ii}));
    this_params = this_params.fit_params;
    h_params(ii,:) = this_params('CI',:);
end

med_h_params = median(h_params{:,:});



%% plot posterior on p_common vs actual for different subjects
% just going to do this for Juno right now


A_sig = mean(med_j_params(2:3));
V_sig = med_j_params(1);
prior_sig = med_j_params(4);
p_common = med_j_params(5);
% make array of predicted post common

pred_post_common = perc_sing_sac(1:20,[1:2 4]);%lil hacky
pred_post_common.pred = zeros(1,20)';
for ii = 1:height(pred_post_common)
    pred_post_common(ii,:).pred = get_post_common(pred_post_common.A_tar(ii),pred_post_common.V_tar(ii),0,A_sig,V_sig,prior_sig,p_common);
end


[g, g_disp]= findgroups(pred_post_common(:,'disp'));
pred_means = splitapply(@mean,pred_post_common.pred,g);
figure(2)
hold on

plot(g_disp.disp,pred_means,'k')
legend('Actual','posterior')

saveas(gcf,'results\behavior\Juno_psing_modeled','png');
%% do again with humans

A_sig = mean(med_h_params(2:3));
V_sig =med_h_params(1);
prior_sig = med_h_params(4);
p_common = med_h_params(5);
% make array of predicted post common
pred_post_common = perc_sing_sac(1:20,[1:2 4]);%lil hacky
pred_post_common.pred = zeros(1,20)';
for ii = 1:height(pred_post_common)
    pred_post_common(ii,:).pred = get_post_common(pred_post_common.A_tar(ii),pred_post_common.V_tar(ii),0,A_sig,V_sig,prior_sig,p_common);
end

%this works much better

[g, g_disp]= findgroups(pred_post_common(:,'disp'));
pred_means = splitapply(@mean,pred_post_common.pred,g);
figure(1)
hold on

plot(g_disp.disp,pred_means,'k')
legend('Actual','posterior')
saveas(gcf,'results\behavior\Human_psing_modeled','png');

