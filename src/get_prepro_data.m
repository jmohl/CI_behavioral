%% preprocess data to get condition and response vectors
%
% -------------------
% Jeff Mohl
% 2/15/19
% -------------------
%
% Description: preprocess the data to work with fitmodel format.
%
% Output format:
% Model [x 1 x]: unity judgement, n saccades (max =2)
% Model [1,2,x]: localization, bayes V1 CxV array with hist counts in each
% bin. Saccades lumped together into one distribution
% Model [2,1,x]: localization, bayes V2 CxVxA array where the counts are
% taken in a grid defined by [Vcoord,Acoord] of a two saccade pair. If one
% saccade Vcoord = Acoord

% TODO: some cleaning/formatting of this code for readability.
function [conditions,responses] = get_prepro_data(data,model,fitoptions)

% preprocess data to work with fitmodel
% get condition vectors
data = data(strcmp(data.trial_type,'AV'),:);
conditions = table2array(unique(data(:,{'A_tar','V_tar'}),'rows'));
midpoints = fitoptions.eval_midpoints;
edges = fitoptions.eval_range;

unity_judge = model(2) == 1;
location_estimate = model(2) == 2;
if model(2) == 3 %do joint fit
    unity_judge = 1;
    location_estimate = 1;
end

if unity_judge %unity judgement model
    %convert into limited table %todo maybe not this?
    data_array = data(:,[2,4,5,16]); %tr num, Atar, Vtar, nsaccades
    data_array = table2array(data_array);
    %currently dealing with more than 1 saccade by saying it is just 2. this
    %might change in the future but I need to carefully look at some of the
    %nsac code in tidy_data project to see what to do about it
    data_array(data_array(:,4)>1,4) = 2;
    
    respbins = unique(data_array(:,4));
    responses_u = zeros(length(conditions),length(respbins));
    for ic = 1:length(responses_u)
        % get reponse for each condition (1 column single counts, 1 double
        % counts)
        responses_u(ic,1) = sum(data_array(:,2) == conditions(ic,1)&data_array(:,3) == conditions(ic,2) & data_array(:,4) == 1); %count single saccade trials
        responses_u(ic,2) = sum(data_array(:,2) == conditions(ic,1)&data_array(:,3) == conditions(ic,2) & data_array(:,4) == 2); %count double saccade trials
    end
end
if location_estimate %localization
    if model(1) == 1 %bayes optimal, V1, all saccades lumped together
        %get responses for target localization, binned in 1 degree bins.
        respbins = fitoptions.eval_range;
        responses_L = zeros(length(conditions),length(midpoints));
        for ic = 1:length(conditions)
            this_data = data(data{:,{'A_tar'}}==conditions(ic,1) & data{:,{'V_tar'}} == conditions(ic,2),:);
            %get only valid saccades for each trial
            valid_sacs = this_data.valid_endpoints;
            valid_sacs = vertcat(valid_sacs{:});
            valid_sacs = valid_sacs(:,1); %only including xcoord
            responses_L(ic,:) = histcounts(valid_sacs',respbins); %count single saccade trials
        end
    else %bayes optimal, v2, split auditory and visual saccades
        responses_L = zeros(length(conditions),length(midpoints),length(midpoints)); %[condition x Vlocs x Alocs]
        for ic = 1:length(conditions)
            this_data = data(data{:,{'A_tar'}}==conditions(ic,1) & data{:,{'V_tar'}} == conditions(ic,2),:);
            %get auditory and visual saccades for each trial
            A_sacs = this_data.A_endpoints;
            V_sacs = this_data.V_endpoints;
            %get single saccade locations
            single_sacs = this_data(this_data.n_sacs == 1,:).valid_endpoints;
            %for all single saccade trials, A_sac = V_sac
            A_sacs(this_data.n_sacs == 1) = single_sacs;
            V_sacs(this_data.n_sacs == 1) = single_sacs;
            %format into vectors
            A_sacs = cell2mat(A_sacs);
            V_sacs = cell2mat(V_sacs);
            responses_L(ic,:,:) = histcounts2(V_sacs(:,1),A_sacs(:,1),edges,edges);
        end
    end
end

if unity_judge && location_estimate
    responses{1} = responses_u;
    responses{2} = responses_L;
elseif unity_judge
    responses = responses_u;
elseif location_estimate
    responses = responses_L;
end

end