%% get response endpoints for a given condition
%
% -------------------
% Jeff Mohl
% 5/31/18
% -------------------
%
% Description: this function will get a vector of saccade endpoints for a
% specified condition. These endpoints are filtered to be only those that
% occur within the specified time window (usually after go time, before end
% time+100 ms buffer)
%
% Valid saccades meet the following description:
% 1: occurs between go_time and end_time + buffer_time (optional buffer time)
% 2: if more than one saccade, saccades are separated by at least 3 degrees
% (if less than 3 degrees, first saccade is considered 'corrected' and only
% second saccade is used)
%
% INPUTS:
%   -data - tidy_data structure for desired trial type.
% OUTPUTS:
%   -endpoints - saccade endpoints [xcoord, ycoord] for every saccade
%   occuring in data within the go_time to end_time + buffer intervals for each trial
%   - A_endpoints - value of endpoints that most closely matches A target
%   - V_endpoints - value of endpoints that most closely matches V target,
%   if A and V endpoints are the same the A_endpoints value is reassigned
%   as the second closest.

function [endpoints, A_endpoints,V_endpoints] = get_response_endpoints(data, req_fix, buffer_time)
%parameters
if nargin < 2
    req_fix = true; %only run on trials where fixation was held until go cue
end

if nargin < 3
    buffer_time = 0; %in ms, amount of time to include after stop time. For saccades initiated before end of trial.
end

if req_fix
    data = data(~isnan(data.go_time),:);
    if isempty(data)
        endpoints=[];
        return
    end
end

endpoints = cell(height(data),1);
A_endpoints = endpoints;
V_endpoints = endpoints;
for i=1:height(data)
    this_endpoints = data.sac_endpoints{i}; %for saccade data
    this_sac_interval = data.sac_intervals{i};
    if ~isempty(this_endpoints)
        interval_inds = this_sac_interval(:,2) > data.go_time(i) & this_sac_interval(:,2) < (data.end_time(i) + buffer_time); %note that endpoints(:,2) contains the endpoint time of the saccade
        this_endpoints = this_endpoints(interval_inds,:);
        %for auditory and visual saccades, only count the closest saccade
        if strcmp(data.trial_type(i),'A') | strcmp(data.trial_type(i),'V')
            tar_loc = max(data(i,:).A_tar,data(i,:).V_tar);
            [~,closest_ind] = min(abs(this_endpoints(:,1)-tar_loc));
            endpoints{i} = this_endpoints(closest_ind,1:2);
        else
            % for multi-target case return all the endpoints, but also
            % split endpoints into A and V distributions by assigning these
            % as the closest saccade to each target.
            
            % BIG THING: I'm going to treat saccades that are less than 3
            % degrees apart as corrective rather than as two response
            % saccades. This is going to have the effect of increaasing the
            % number of single saccade trials, but I think that is more
            % representative of intent. 
            
            if size(this_endpoints,1) > 1
                corrective_ind = zeros(size(this_endpoints,1),1);
                for jj = 1:size(this_endpoints,1)-1
                    corrective_ind(jj) = abs(this_endpoints(jj,1) - this_endpoints(jj+1,1)) <= 3;
                end
                this_endpoints = this_endpoints(~corrective_ind,:);
            end
            endpoints{i} = this_endpoints;
           
            if size(this_endpoints,1) > 1
                A_tar = data.A_tar(i);V_tar = data.V_tar(i);
                [~,A_closest_ind] = sort(abs(this_endpoints(:,1)-A_tar));
                A_endpoints{i} = this_endpoints(A_closest_ind(1),1:2);
                [~,V_closest_ind] = min(abs(this_endpoints(:,1)-V_tar));
                V_endpoints{i} = this_endpoints(V_closest_ind,1:2);
                if A_endpoints{i} == V_endpoints{i} %V is assumed to be more accurate, so if there is one saccade that is "closest" to both targets, A is changed to the 2nd closest.
                    A_endpoints{i} = this_endpoints(A_closest_ind(2),1:2);
                end
            end
        end
    else
        endpoints{i} = [];
    end
end

end
