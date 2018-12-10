%% find percentage of single saccade trials per condition
%
% -------------------
% Jeff Mohl
% 12/10/18
% -------------------
%
% Description: for each condition, determine the fraction of trials which
% have a single saccade (indicating a single target was perceived). This
% uses the n_sacs field in tidy data, which counts saccades which occur
% between the go cue and the end state (plus a buffer of around 100ms to
% account for saccades which were initiated before the end state).
%
%


function results = get_perc_sing_sac(data)
AV_data = data(strcmp(data.trial_type, 'AV'),:);
sac_array = AV_data(:,{'A_tar','V_tar','n_sacs'});
sac_array.disp = sac_array.A_tar - sac_array.V_tar;
[g, pairs] = findgroups(sac_array(:,1:2));
p_single = [];
for i=1:20
    p_single(i) = sum(sac_array.n_sacs(g==i) == 1)/length(sac_array.n_sacs(g==i));%percent single saccade
end
results = pairs;
results.psingle = p_single';
results.disp = results.A_tar - results.V_tar;
results.ID = data.file_ID(1:20);

end