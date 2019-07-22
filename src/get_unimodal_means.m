%% estimate mean values from unimodal trials, grouped by target
%
% -------------------
% Jeff Mohl
% 7/22/19
% -------------------
%
% Description: This function finds the mean values for the unimodal
% localization trials (split by target value) and returns a 'condition'
% vector of the same size and format. This can be fed into the model to
% make the fits use the mean responses to a given target location instead
% of the actual target locations for purposes of modeling.

function [mean_resps] = get_unimodal_means(conditions, data,model)

% get A trials and V trials
A_data = data(strcmp(data.trial_type,'A'),:);
V_data = data(strcmp(data.trial_type,'V'),:);
[gA,vA] = findgroups(A_data.A_tar);
[gV,vV] = findgroups(V_data.V_tar);

A_sacs = vertcat(A_data.valid_endpoints{:});
V_sacs = vertcat(V_data.valid_endpoints{:});
A_means = splitapply(@mean, A_sacs,gA);
V_means = splitapply(@mean, V_sacs,gV);
A_means = A_means(:,1); %only care about x coord
V_means = V_means(:,1);

if model(3) == 4 %unisensory localization model
    mean_resps{1} = A_means;
    mean_resps{2} = V_means;
else %all the multi-sensory, multi-target models
    mean_resps = conditions;
    for ind = 1:length(conditions)
        mean_resps(ind,1) = A_means(vA == conditions(ind,1));
        mean_resps(ind,2) = V_means(vV == conditions(ind,2));
    end
end

end