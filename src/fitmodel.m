%% Fit model to data
%
% -------------------
% Jeff Mohl
% 2/1/19
% -------------------
%
% Description: this is code for performing the model fitting steps for the
% given model and given subject. Will return the best fit parameters, the
% nll for the model fit... maybe more in the future
%
% Inputs:
% data(1) trial number
% data(2) A tar
% data(3) V tar
% data(4) response (number of saccades for unity judgement case)
% theta - contains model fit parameters  
% model(1) Model type (1=Bayesian, 2=null, 3=switching)
% model(2) Response type (1=unity judgement, 2 = location)
% model(3) estimation proceedure (1=numerical integration)

%% hack code for starting this process
global MAXRNG
MAXRNG = 60;

cd('C:\Users\jtm47\Documents\Projects\CI_behavioral')
addpath('data','src','src\lautils', 'src\plotting');

model = [1 2 1];

data = load('H08_AVD2_2018_08_10_tidy.mat');
data= data.tidy_data;
%get only AV trials
data = data(strcmp(data.trial_type,'AV') & ~isnan(data.go_time),:);
data = data(data.n_sacs > 0,:);%only include trials with a saccade

%setting test parameters
theta = [5,5,5,.5,.1]; %v_sig, A_sig, prior_sig, prior_common, lambda (lapse rate)

%setting fitting procedure options
fmin_options = optimset('MaxFunEvals',1000,'MaxIter',500,'Display','iter','TolX',1e-3);

eval_n = MAXRNG*2+1;
eval_range = linspace(-MAXRNG,MAXRNG,eval_n);
eval_midpoints = linspace(-60+60/eval_n,60-60/eval_n,eval_n-1);


debug = 1;
%% process data for fitting procedure
% get condition vectors
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
    respbins = eval_range;
    responses = zeros(length(conditions),length(eval_midpoints));
    for ic = 1:length(conditions)
        this_data = data(data{:,{'A_tar'}}==conditions(ic,1) & data{:,{'V_tar'}} == conditions(ic,2),:);
        %get only valid saccades for each trial
        valid_sacs = get_response_endpoints(this_data,1,100);
        valid_sacs = vertcat(valid_sacs{:});
        valid_sacs = valid_sacs(:,1); %only including xcoord
        responses(ic,:) = histcounts(valid_sacs',respbins); %count single saccade trials
    end
end


%% start of actual fitting procedure
%adapting from Acerbi, this fitting procedure will progress in 2 steps to
%start. First, a grid with 1000 points spanning reasonable parameter space will be
%used to determine starting points for the model. The best 5 starting
%points will be used as initial points for fminsearch optimization, with
%parameter and nll values saved out.

%todos: validation, multi-step fitting

%2/1/19 note: used timeit to test this likelihood function as 0.0091 sec
%for human subject H08 unity judgement. Unfortunately for the localization
%judgement case it takes 2.28 seconds (with eval_n = 250, with eval_n = 100
%significantly faster at 0.147 sec), which is pretty slow considering the
%number of times this function needs to be called. Also the likelihood does
%not converge over 500 iterations, which takes 2 minutes.
%
% f=@() datalike_minsearch(theta);
% timeit(f)

datalike_minsearch = @(theta)datalike(conditions,responses,theta,model,eval_midpoints);
tic
[fit_results,fit_results_nll,~,~] = fminsearch(datalike_minsearch,theta,fmin_options);
toc
%this was able to run the optimization procedure in 6.4 sec for the unity
%judgement case. For the localization case it takes much longer.




if debug 
    %plot some things for comparing with behavior
    if model(2) == 1
    [nll,modelfit] = datalike(conditions,responses,fit_results,model,eval_midpoints);
    plot_psingle(data,modelfit);
    end
    if model(2) == 2
        [nll,modelfit] = datalike(conditions,responses,fit_results,model,eval_midpoints);
        for ic = 1:length(conditions)
            %plotting the real saccade distributions and those predicted by
            %the model for every condition
            plot_modelhist(responses(ic,:),modelfit(ic,:),eval_midpoints)
            title(sprintf('%d A %d V',conditions(ic,:)))
            set(gcf,'Position',[100,60,1049,895])
        end
    end
        
end










