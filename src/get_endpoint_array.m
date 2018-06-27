%% get endpoints and target locs in cell array form for fitting
%
% -------------------
% Jeff Mohl
% 6/27/18
% -------------------
%
% Description: this repackages saccade endpoint data for a given trial type
% in the cell array format needed for efficient model fitting.

function [endpoint_data,tars] = get_endpoint_array(data,type, req_fix,sac_buffer)

% create endpoint matrix
% only using data from AV trial type
if strcmp(type,'AV')
    this_data = data(strcmp(data.trial_type,'AV'),:);
    tar_pairs = sortrows(unique([this_data.A_tar,this_data.V_tar],'rows')); %first column is A_tar, second is V_tar
else
    this_data = data(strcmp(data.trial_type,type),:);
    %first column is A_tar, second is V_tar
    if strcmp(type,'A')
        tar_pairs(:,1) = sortrows(unique(this_data.A_tar(~isnan(this_data.A_tar))));
        tar_pairs(:,2) = NaN;
    else
        tar_pairs(:,2) = sortrows(unique(this_data.V_tar(~isnan(this_data.V_tar))));
        tar_pairs(:,1) = NaN;
    end
    
end
% get AV endpoints for all target pairs, and save in cell array for faster
% fitting
for i=1:length(tar_pairs)
    A_tar = tar_pairs(i,1);
    V_tar = tar_pairs(i,2);
    if strcmp(type,'AV')
        endpoints = get_response_endpoints(this_data(this_data.A_tar == A_tar & this_data.V_tar == V_tar,:),req_fix,sac_buffer); % require fix, 100ms buffer, checked working 6/22
    else
        endpoints = get_response_endpoints(this_data(this_data.A_tar == A_tar | this_data.V_tar == V_tar,:),req_fix,sac_buffer); % require fix, 100ms buffer, checked working 6/22
    end
    %trim to only horizontal component
    endpoint_data{i} = endpoints(:,1);
end
tars = tar_pairs(:,all(~isnan(tar_pairs))); %trim nans for unimodal cases.
end