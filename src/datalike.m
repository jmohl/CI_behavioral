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
% data(2) A tar
% data(3) V tar
% data(4) response (number of saccades for unity judgement case)
% theta - contains model fit parameters  
% model (?) todo
%

%% start of actual likelihood fitting code
function [nll,prmat_unity] = datalike(data,theta)

global MAXRNG
debug = 0;

V_sig = theta(1);%.V_sig; %vis target sigma
A_sig = theta(2);%.A_sig; %close aud target sigma
prior_sig = theta(3);%.prior_sig;%sigma of centrality prior
p_common = theta(4);%.prior_common;%prior on common cause
lambda = theta(5); %lapse probability

%HACK fminsearch does not allow bounds, so I'm including this really hacky
%way to require that the likelihood is not calculated for impossible
%values. This will not be necessary when I switch to using either bads or
%fminsearchbnd, both of which are third party.
if min(theta(1:3)) <= 0.5 || p_common > 1 || p_common < 0 
    nll = 10000;
    return;
end

xrange = linspace(-60,60,250); %range for calculating likelihoods.
xrange_V(1,:,1) = xrange;
xrange_A(1,1,:) = xrange';

% get condition vectors
conds = unique(data(:,2:3),'rows');
conds_A = conds(:,1);
conds_V = conds(:,2);
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
likec2 = bsxfun(@times, int_V, int_A);

% calculate the C=1 PDF
%NOTE: this code is adapted from lines 399:410 in
%VestBMS_BimodalLeftRightDatalike.m by luigi acerbi
mutilde = bsxfun(@plus, xrange_A.*A_sig^2, xrange_V.*V_sig^2)./(A_sig^2 + V_sig^2);
sigma2tilde = A_sig^2.*V_sig^2./(A_sig^2 + V_sig^2);
mucdf = mutilde.*prior_sig^2./(sigma2tilde + prior_sig^2);
sigmacdf = sqrt(sigma2tilde./(sigma2tilde + prior_sig^2))*prior_sig;
intc1 = (bsxfun_normcdf(MAXRNG, mucdf, sigmacdf) - bsxfun_normcdf(-MAXRNG, mucdf, sigmacdf));
likec1 = intc1 .* bsxfun_normpdf(xrange_A,xrange_V,sqrt(A_sig^2 + V_sig^2)) .* ...
    bsxfun_normpdf(mutilde,0,sqrt(sigma2tilde +prior_sig^2));

% get posterior for (C=1|xv,xa)
c1post = (likec1 * p_common)./(likec1 * p_common + (1-p_common)*likec2);

%judgement rule, if post > 0.5, choose unity
w1_unity = zeros(size(c1post));
w1_unity(c1post > 0.5) = 1;
w1_unity(c1post == 0.5) = 0.5;

%plots for debugging
if debug
    figure;imagesc(squeeze(c1post));
    figure;imagesc(squeeze(w1_unity));
end

%integrate over xa xv values for every condition in condition vector
% try both using qtrapz and VestBMS_finalqtrapz. Again adapted from acerbi
xpdf_V = bsxfun_normpdf(xrange_V, conds_V,V_sig);
xpdf_V = bsxfun(@rdivide, xpdf_V, qtrapz(xpdf_V, 2)); % Not multiplying by volume element

xpdf_A = bsxfun_normpdf(xrange_A, conds_A, A_sig);
xpdf_A = bsxfun(@rdivide, xpdf_A, qtrapz(xpdf_A, 3));  % Not multiplying by volume element(JM?)
prmat_unity = zeros(numel(conds_V), 2);
%prmat_unity(:,1) = VestBMS_finalqtrapz(xpdf_V,xpdf_A,w1_unity);    % Not multiplying by volume element (xpdfs did not)
prmat_unity(:,1) = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), w1_unity), 2), 3); %this is a recreation of what VestBMS_finalqtrapz does but not done in C. 
prmat_unity(:,2) = 1 - prmat_unity(:,1);

% Fix probabilities
%JM p cant be greater than 1 or smaller than 0.
prmat_unity = min(max(prmat_unity,0),1);

%add in change for random response
prmat_unity = lambda/2 + (1-lambda)*prmat_unity;
if debug
   %plot decision rule by target sep, for rough comparison with data (will
   %actually depend on eccentricity as well (which is why I'm plotting
   %dots, there will be a few dots for each sep value
   tar_sep = abs(conds_V - conds_A);
   figure;plot(tar_sep,prmat_unity(:,1),'k.')
   g=findgroups(tar_sep);
   sep_means = splitapply(@mean,prmat_unity(:,1),g);
   hold on; plot(unique(tar_sep),sep_means);
   legend('individual pairs','mean by sep')
   title('p single by sep, predicted by model')
   xlabel('degrees of separation')
   ylabel('p single saccade response')
end

% calculate negative log likelihood of data
nll = -1*sum(responses.*log(prmat_unity));
nll = sum(nll); %sum across bins


end
