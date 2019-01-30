%% Likelihood calculation for various models
%
% -------------------
% Jeff Mohl
% 1/30/19
% -------------------
%
% Description: this the updated code for calculating the log likelihood of
% a given dataset under each model. Initially making this work for the
% unity judgement case but will hopefully expand to work with all.
%
% Inputs:
% data(1) trial number
% data(2) task modality 
% model name of model being used(?)
%
%

%% hack code for starting this process

data = load('H08_AVD2_2018_08_10_tidy.mat');
data= data.tidy_data;
%get only AV trials
data = data(strcmp(data.trial_type,'AV') & ~isnan(data.go_time),:);
%convert into limited table
data = data(:,[2,4,5,16]); %tr num, Atar, Vtar, nsaccades
data = table2array(data);
%currently dealing with more than 1 saccade by saying it is just 2. this
%might change in the future but I need to carefully look at some of the
%nsac code in tidy_data project to see what to do about it 
data(data(:,4)>1,4) = 2;

%setting test parameters
theta = [1,3,5,.5,.1]; %v_sig, A_sig, prior_sig, prior_common, lambda (lapse rate)

%% start of actual likelihood fitting code

V_sig = theta(1);%.V_sig; %vis target sigma
A_sig = theta(2);%.A_sig; %close aud target sigma
prior_sig = theta(3);%.prior_sig;%sigma of centrality prior
p_common = theta(4);%.prior_common;%prior on common cause
lamba = theta(5); %lapse probability

xrange = -45:.5:45; %range for calculating likelihoods.

% get condition vectors

% get reponse vectorx2 for each condition (1 column single counts, 1 double
% counts )

% calculate the C=2 PDF at each value of xrange, joint probability of xa xv

% calculate the C=1 PDF

% get posterior for (C=1|xv,xa)
c1post = (c1like * p_common)/(c1like * p_common + (1-p_common)*c2like);

%integrate over xa xv values for every condition in condition vector
% try both using qtrapz and VestBMS_finalqtrapz

% Fix probabilities
prmat = min(max(prmat,0),1); %JM p cant be greater than 1 or smaller than 0.
prmat_unity = min(max(prmat_unity,0),1);

% calculate log likelihood of data
ll = sum(resp.*log(prmat_unity));
