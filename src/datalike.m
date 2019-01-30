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
theta = [5,5,5,.5,.1]; %v_sig, A_sig, prior_sig, prior_common, lambda (lapse rate)
MAXRNG = 45;
realmin = 0.00001;

%% start of actual likelihood fitting code

V_sig = theta(1);%.V_sig; %vis target sigma
A_sig = theta(2);%.A_sig; %close aud target sigma
prior_sig = theta(3);%.prior_sig;%sigma of centrality prior
p_common = theta(4);%.prior_common;%prior on common cause
lamba = theta(5); %lapse probability

xrange = -60:.5:60; %range for calculating likelihoods.
xrange_V = xrange;
xrange_A = xrange';

% get condition vectors
conds = unique(data(:,2:3),'rows');
respbins = unique(data(:,4));

% get reponse for each condition (1 column single counts, 1 double
% counts )
responses = zeros(length(conds),length(respbins));
for ic = 1:length(responses)
    responses(ic,1) = sum(data(:,2) == conds(ic,1)&data(:,3) == conds(ic,2) & data(:,4) == 1); %count single saccade trials in cond ci
    responses(ic,2) = sum(data(:,2) == conds(ic,1)&data(:,3) == conds(ic,2) & data(:,4) == 2); %count double saccade trials in cond ci
end

% calculate the C=2 PDF at each value of xrange, joint probability of xa xv
%NOTE: this code is adapted from lines 357:371 in
%VestBMS_BimodalLeftRightDatalike.m by luigi acerbi
%NOTE2: bsxfun_normpdf and _normcdf are also written by acerbi, but I have
%copied them into this project src
muc2_V = xrange_V.*prior_sig^2/(V_sig^2 + prior_sig^2);
sigmac2_V = V_sig*prior_sig/sqrt(V_sig^2 + prior_sig^2);
muc2_A = xrange_A.*prior_sig^2/(A_sig^2 + prior_sig^2);
sigmac2_A = A_sig*prior_sig/sqrt(A_sig^2 + prior_sig^2);
int_V = (bsxfun_normcdf(MAXRNG,muc2_V,sigmac2_V) - bsxfun_normcdf(-MAXRNG,muc2_V,sigmac2_V));
int_A = (bsxfun_normcdf(MAXRNG,muc2_A,sigmac2_A) - bsxfun_normcdf(-MAXRNG,muc2_A,sigmac2_A));
int_V = int_V .* bsxfun_normpdf(xrange_V,0,sqrt(V_sig^2 + prior_sig^2));
int_A = int_A .* bsxfun_normpdf(xrange_A,0,sqrt(A_sig^2 + prior_sig^2));
likec2 = bsxfun(@times, int_V, int_A) + realmin;

% calculate the C=1 PDF
%NOTE: this code is adapted from lines 399:410 in
%VestBMS_BimodalLeftRightDatalike.m by luigi acerbi
mutilde = bsxfun(@plus, xrange_A.*A_sig^2, xrange_V.*V_sig^2)./(A_sig^2 + V_sig^2);
sigma2tilde = A_sig^2.*V_sig^2./(A_sig^2 + V_sig^2);
mucdf = mutilde.*prior_sig^2./(sigma2tilde + prior_sig^2);
sigmacdf = sqrt(sigma2tilde./(sigma2tilde + prior_sig^2))*prior_sig;
intc1 = (bsxfun_normcdf(MAXRNG, mucdf, sigmacdf) - bsxfun_normcdf(-MAXRNG, mucdf, sigmacdf));
likec1 = intc1 .* bsxfun_normpdf(xrange_A,xrange_V,sqrt(A_sig^2 + V_sig^2)) .* ...
    bsxfun_normpdf(mutilde,0,sqrt(sigma2tilde +prior_sig^2)) + realmin;

% get posterior for (C=1|xv,xa)
c1post = (c1like * p_common)/(c1like * p_common + (1-p_common)*c2like);

%integrate over xa xv values for every condition in condition vector
% try both using qtrapz and VestBMS_finalqtrapz

% Fix probabilities
prmat = min(max(prmat,0),1); %JM p cant be greater than 1 or smaller than 0.
prmat_unity = min(max(prmat_unity,0),1);

% calculate log likelihood of data
ll = sum(resp.*log(prmat_unity));
