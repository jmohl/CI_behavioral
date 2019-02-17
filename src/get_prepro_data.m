%% preprocess data to get condition and response vectors
%
% -------------------
% Jeff Mohl
% 2/15/19
% -------------------
%
% Description: preprocess the data to work with fitmodel format.
%

function [conditions,responses] = get_prepro_data(data,model,fitoptions)

% preprocess data to work with fitmodel
% get condition vectors
data = data(strcmp(data.trial_type,'AV'),:);
conditions = table2array(unique(data(:,{'A_tar','V_tar'}),'rows'));

if model(2) == 1
    %convert into limited table %todo maybe not this?
    data = data(:,[2,4,5,16]); %tr num, Atar, Vtar, nsaccades
    data = table2array(data);
    %currently dealing with more than 1 saccade by saying it is just 2. this
    %might change in the future but I need to carefully look at some of the
    %nsac code in tidy_data project to see what to do about it
    data(data(:,4)>1,4) = 2;
    
    respbins = unique(data(:,4));
    responses = zeros(length(conditions),length(respbins));
    for ic = 1:length(responses)
        % get reponse for each condition (1 column single counts, 1 double
        % counts)
        responses(ic,1) = sum(data(:,2) == conditions(ic,1)&data(:,3) == conditions(ic,2) & data(:,4) == 1); %count single saccade trials
        responses(ic,2) = sum(data(:,2) == conditions(ic,1)&data(:,3) == conditions(ic,2) & data(:,4) == 2); %count double saccade trials
    end
else
    %get responses for target localization, binned in 1 degree bins.
    respbins = fitoptions.eval_range;
    responses = zeros(length(conditions),length(fitoptions.eval_midpoints));
    for ic = 1:length(conditions)
        this_data = data(data{:,{'A_tar'}}==conditions(ic,1) & data{:,{'V_tar'}} == conditions(ic,2),:);
        %get only valid saccades for each trial
        valid_sacs = this_data.valid_endpoints;
        valid_sacs = vertcat(valid_sacs{:});
        valid_sacs = valid_sacs(:,1); %only including xcoord
        responses(ic,:) = histcounts(valid_sacs',respbins); %count single saccade trials
    end
end
end