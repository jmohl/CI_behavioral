%% preprocess data to get condition and response vectors
%
% -------------------
% Jeff Mohl
% 2/15/19
% -------------------
%
% Description: preprocess the data to work with fitmodel format.
%
% Output format is different depending on model, for each condition (C): 
% Model [x x 1]: unity judgement task, Cx2 array of number of one saccade
% and two saccade trials per condition 
% Model [1,x,2]: localization only, Cx100x100 array of saccade endpoints 
% where the counts are taken in a grid defined by [Vcoord,Acoord] of a two
% saccade pair. If one saccade Vcoord = Acoord Model
% Model [1,x,3]: Joint fit, cell array containing both of the above
% response matrices.

% TODO: update for dynamic binning, return A edges and V edges, as well as
% midpoints.
function [conditions,responses] = get_prepro_data(data,model)
global fitoptions MAXRNG
% preprocess data to work with fitmodel

switch model(3)
    case 1
        unity_judge = 1;
        location_estimate = 0;
        unisensory_loc = 0;
        
    case 2
        unity_judge = 0;
        location_estimate = 1;
        unisensory_loc = 0;
        
    case 3
        unity_judge = 1;
        location_estimate = 1;
        unisensory_loc = 0;
        
    case 4
        unity_judge = 0;
        location_estimate = 0;
        unisensory_loc = 1;
end

% get condition vectors
if ~unisensory_loc
data = data(strcmp(data.trial_type,'AV'),:);
conditions = table2array(unique(data(:,{'A_tar','V_tar'}),'rows'));
else
    A_data = data(strcmp(data.trial_type,'A'),:);
    V_data = data(strcmp(data.trial_type,'V'),:);
    A_conditions = table2array(unique(A_data(:,{'A_tar'}),'rows'));
    V_conditions = table2array(unique(V_data(:,{'V_tar'}),'rows'));
end
midpoints = fitoptions.eval_midpoints;

if unity_judge %unity judgement model
    %convert into limited table %todo maybe not this?
    data_array = data(:,{'trial','A_tar','V_tar','n_sacs'}); %tr num, Atar, Vtar, nsaccades
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
    if fitoptions.strict_filter
        data = strict_single_filter(data);
        data = data(logical(data.valid_tr),:);
    end
    
    if fitoptions.dynamic_bins
        responses_L = zeros(length(conditions),28,28); %[condition x Vlocs x Alocs]
    else
        responses_L = zeros(length(conditions),length(midpoints),length(midpoints));
        Aedges = fitoptions.eval_range;
        Vedges = fitoptions.eval_range;%range in degrees for evaluation
    end
    for ic = 1:length(conditions)
        V_tar = conditions(ic,2);
        A_tar = conditions(ic,1);
        this_data = data(data{:,{'A_tar'}}==A_tar & data{:,{'V_tar'}} == V_tar,:);
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
        if fitoptions.dynamic_bins %working on
            Vedges = [-MAXRNG:6:V_tar - 12, V_tar - 9:3:V_tar - 6, V_tar-5:1:V_tar + 5, V_tar + 6:3:V_tar + 9, V_tar + 12:6:MAXRNG]; %probably a nicer way to write this.
            Aedges = [-MAXRNG:6:A_tar - 12, A_tar - 9:3:A_tar - 6, A_tar-5:1:A_tar + 5, A_tar + 6:3:A_tar + 9, A_tar + 12:6:MAXRNG]; %probably a nicer way to write this.
        end
%         fitoptions.eval_A = Aedges;
%         fitoptions.eval_V = Vedges;
        responses_L(ic,:,:) = histcounts2(V_sacs(:,1),A_sacs(:,1),Vedges,Aedges);
%          %for dynamic bin testing purposes
%         figure
%         hist3([V_sacs(:,1),A_sacs(:,1)],'Edges',{Vedges,Aedges});
%         xlabel('V saccades')
%         ylabel('A saccades')
%         title(sprintf('%d A %d V',A_tar,V_tar))
    end
end

% unisensory localization
if unisensory_loc
        responses_A = zeros(length(A_conditions),length(midpoints));
        responses_V = zeros(length(V_conditions),length(midpoints));
        Aedges = fitoptions.eval_range;
        Vedges = fitoptions.eval_range;%range in degrees for evaluation
        
    for ia = 1:length(A_conditions)
        A_tar = A_conditions(ia);
        this_data = A_data(A_data{:,{'A_tar'}}==A_tar,:);
        %get auditory  saccades for each trial
        A_sacs = this_data.valid_endpoints;
        %format into vectors
        A_sacs = cell2mat(A_sacs);
        responses_A(ia,:) = histcounts(A_sacs(:,1),Aedges);
    end
    for iv = 1:length(V_conditions)
        V_tar = V_conditions(iv);
        this_data = V_data(V_data{:,{'V_tar'}}==V_tar,:);
        %get auditory  saccades for each trial
        V_sacs = this_data.valid_endpoints;
        %format into vectors
        V_sacs = cell2mat(V_sacs);
        responses_V(iv,:) = histcounts(V_sacs(:,1),Vedges);
    end
end

if unisensory_loc
    responses{1} = responses_A;
    responses{2} = responses_V;
    conditions{1} = A_conditions;
    conditions{2} = V_conditions;
elseif unity_judge && location_estimate
    responses{1} = responses_u;
    responses{2} = responses_L;
elseif unity_judge
    responses = responses_u;
elseif location_estimate
    responses = responses_L;
end

end