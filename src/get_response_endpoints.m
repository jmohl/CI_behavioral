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
% INPUTS:
%   -data - tidy_data structure for desired trial type.
% OUTPUTS: 
%   -endpoints - saccade endpoints [xcoord, ycoord, time] for every saccade
%   occuring in data within the go_time to end_time + buffer intervals for each trial


function endpoints = get_response_endpoints(data, req_fix, buffer_time)
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

for i=1:height(data)
    this_endpoints = data.sac_endpoints{i}; %for saccade data
    if ~isempty(this_endpoints)
        interval_inds = this_endpoints(:,3) > data.go_time(i) & this_endpoints(:,3) < (data.end_time(i) + buffer_time); %note that endpoints(:,3) contains the time point for the end of the saccade
        %for auditory and visual saccades, only count first saccade
        if strcmp(data.trial_type(i),'A') | strcmp(data.trial_type(i),'V')
            endpoints{i} = this_endpoints(find(interval_inds,1,'first'),1:3);
        else
            endpoints{i} = this_endpoints(interval_inds,1:3);
        end
    else
        endpoints{i} = [];
    end
end

end
