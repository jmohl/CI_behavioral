%% calculate stats for a given dataset
%
% -------------------
% Jeff Mohl
% 7/24/19
% -------------------
%
% Description: generates plots of raw behavioral data for a few demo
% conditions. The purpose of these plots is to show that people and monkeys
% basically perform the task. There are 3 plots, which can be repeated for
% each subject or across subjects. These are:
% 1: unisensory localization accuracy for each target
% 2: saccade traces, with each frame plotted to show speed and fixation etc
% 3: saccade histograms for the same condition, maybe one each for
% integrated/segregated

function raw_behavior_stats = get_raw_behavior_stats(this_data,subject)
    
%only use valid trials
    this_data = this_data(logical(this_data.valid_tr),:); %only include valid trials

    A_data = this_data(strcmp(this_data.trial_type,'A'),:);
    V_data = this_data(strcmp(this_data.trial_type,'V') & abs(this_data.V_tar) ~= 30,:);
    
    [gA,glabA] = findgroups(A_data.A_tar);
    [gV,glabV] = findgroups(V_data.V_tar);
    A_sac = vertcat(A_data.valid_endpoints{:});
    A_sac = A_sac(:,1);
    V_sac = vertcat(V_data.valid_endpoints{:});
    V_sac = V_sac(:,1);
    %get same data from multiple AV conditions
    %ex_conditions 
    AV_ex_data = this_data(strcmp(this_data.trial_type,'AV'),:);
    AV_ex_data = AV_ex_data(abs(AV_ex_data.A_tar - AV_ex_data.V_tar) > 17,:);
    %only want dual saccade data for this stat
    AV_ex_data = AV_ex_data(AV_ex_data.n_sacs > 1,:);
    AV_A_sac = vertcat(AV_ex_data.A_endpoints{:});
    AV_A_sac = AV_A_sac(:,1);
    AV_V_sac = vertcat(AV_ex_data.V_endpoints{:});
    AV_V_sac = AV_V_sac(:,1);
        
    [gAV_A,glabAV_A] = findgroups(AV_ex_data.A_tar);
    [gAV_V,glabAV_V] = findgroups(AV_ex_data.V_tar);
    
    raw_behavior_stats.subject = subject;
    raw_behavior_stats.A_tars = glabA;
    raw_behavior_stats.V_tars = glabV;
    raw_behavior_stats.mean_A = splitapply(@mean,A_sac,gA);
    raw_behavior_stats.mean_V = splitapply(@mean,V_sac,gV);
    raw_behavior_stats.std_A = splitapply(@std,A_sac,gA);
    raw_behavior_stats.std_V = splitapply(@std,V_sac,gV);
    raw_behavior_stats.mean_AV_A = splitapply(@mean,AV_A_sac,gAV_A);
    raw_behavior_stats.mean_AV_V = splitapply(@mean,AV_V_sac,gAV_V);
    raw_behavior_stats.std_AV_A = splitapply(@std,AV_A_sac,gAV_A);
    raw_behavior_stats.std_AV_V = splitapply(@std,AV_V_sac,gAV_V);
    raw_behavior_stats.AV_A_tars = glabAV_A;
    raw_behavior_stats.AV_V_tars = glabAV_V;
    
    %get error under each of these conditions
    raw_behavior_stats.A_error = abs(raw_behavior_stats.mean_A - raw_behavior_stats.A_tars);
    raw_behavior_stats.V_error = abs(raw_behavior_stats.mean_V - raw_behavior_stats.V_tars);
    raw_behavior_stats.AV_A_error = abs(raw_behavior_stats.mean_AV_A - raw_behavior_stats.AV_A_tars);
    raw_behavior_stats.AV_V_error = abs(raw_behavior_stats.mean_AV_V - raw_behavior_stats.AV_V_tars);
    
end