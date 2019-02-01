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
% model (?) todo

global MAXRNG
MAXRNG = 60;

%% hack code for starting this process
cd('C:\Users\jtm47\Documents\Projects\CI_behavioral')
addpath('data','src','src\lautils', 'src\plotting');

data = load('Yoko_AVD2_2019_01_31_tidy.mat');
data= data.tidy_data;
%get only AV trials
data = data(strcmp(data.trial_type,'AV') & ~isnan(data.go_time),:);
data = data(data.n_sacs > 0,:);%only include trials with a saccade
%get valid saccades only, not sure if necessary
%data.valid_endpoints = get_response_endpoints(data,0,100)';

%convert into limited table
data = data(:,[2,4,5,16]); %tr num, Atar, Vtar, nsaccades
data = table2array(data);

%currently dealing with more than 1 saccade by saying it is just 2. this
%might change in the future but I need to carefully look at some of the
%nsac code in tidy_data project to see what to do about it 
data(data(:,4)>1,4) = 2;

%setting test parameters
theta = [5,5,5,.5,.1]; %v_sig, A_sig, prior_sig, prior_common, lambda (lapse rate)

%setting fitting procedure options
fmin_options = optimset('MaxFunEvals',20000,'MaxIter',40000);

debug = 1;

%% start of actual fitting procedure
%adapting from Acerbi, this fitting procedure will progress in 2 steps to
%start. First, a grid with 1000 points spanning reasonable parameter space will be
%used to determine starting points for the model. The best 5 starting
%points will be used as initial points for fminsearch optimization, with
%parameter and nll values saved out.

%todos: validation

%2/1/19 note: used timeit to test this likelihood function as 0.0091 sec
%for human subject H08.
datalike_minsearch = @(theta)datalike(data,theta);
tic
[fit_results,fit_results_nll,~,~] = fminsearch(datalike_minsearch,theta,fmin_options);
toc
%this was able to run the optimization procedure in 6.4 sec

if debug 
    %plot some things for comparing with behavior
    [nll,modelfit] = datalike(data,fit_results);
    plot_psingle(data,modelfit);
end










