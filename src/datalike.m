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
% theta(1) V_sig: variance of visual localization
% theta(2) A_sig: variance of auditory localization
% theta(3) prior_sig: sigma on target location prior
% theta(4) p_common: prior on common cause
% theta(5) lambda: lapse rate, probability of random response occuring
% model(1) Causal inference type (1=bayesian,2=probabilistic fusion, 3 fixed criterion (todo))
% model(2) Response type (1 = unity judgement, 2 = location, 3 = joint)
% model(3) combination rule(1=posterior reweighting, 2 = model selection, 3 = probabilistic fusion (todo, make fusion rate free param?) 4= probability matching (todo))


%% start of likelihood code - currently only working for unity judgement
function [nll,prmat] = datalike(conditions,responses,theta,model,eval_midpoints)

CI_type = model(1);
unity_judge = model(2) == 1;
location_estimate = model(2) == 2;
if model(2) == 3 %do joint fit
    unity_judge = 1;
    location_estimate = 1;
end

combination_rule = model(3); %rule for combining sensory inputs, based on causal judgement 

V_sig = theta(1);%vis target sigma
A_sig = theta(2);%close aud target sigma
prior_sig = theta(3);%sigma of centrality prior
p_common = theta(4);%prior on common cause
lambda = theta(5); %lapse probability
prior_mu = 0; %fixed for now

conds_A = conditions(:,1);
conds_V = conditions(:,2);

%HACK fminsearch does not allow bounds, so I'm including this really hacky
%way to require that the likelihood is not calculated for impossible
%values. In the future I might switch to bads or fminsearchbnd, both of
%which are third party but allow bounds. 
if min(theta(1:3)) <= 0.1 || max(theta(4:5)) > 1 || min(theta(4:5)) < 0
    nll = 1e10;
    return;
end

xrange = eval_midpoints; %range for integration.
xrange_V(1,:,1) = xrange;
xrange_A(1,1,:) = xrange';

%% find the posterior distribution for C = 1 case for all values of xa and xv;

switch CI_type 
    case 1 % bayesian causal inference posterior
        c1post = get_c1post(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,p_common);
    case 2
        c1post = repmat(p_common,1,length(xrange),length(xrange));
end

if unity_judge
    %judgement rule, if post > 0.5, choose unity
    if CI_type == 2 %if probabilistic fusion rule
        w1_unity = c1post;
    else %use fixed decision rule for now
    w1_unity = zeros(size(c1post));
    w1_unity(c1post > 0.5) = 1;
    w1_unity(c1post == 0.5) = 0.5;
    end
end

%% get likelihood functions for A location and V location
if location_estimate
    int_pdf = get_integrate_pdf(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,xrange);
    
    [~,A_seg_pdf,V_seg_pdf] = get_segregate_pdf(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,xrange);
    
    % combine likelihoods weighted by posterior on common cause
    switch combination_rule
        case 1 %if bayesian reweighting
            V_c1c2_pdf = bsxfun(@times,c1post, int_pdf) + bsxfun(@times,(1-c1post), V_seg_pdf);
            A_c1c2_pdf = bsxfun(@times,c1post, int_pdf) + bsxfun(@times,(1-c1post), A_seg_pdf);
        case 2 %if model selection, choose the best model and use that one exclusively (set weight on alternative to 0)
            w1_unity = zeros(size(c1post));
            w1_unity(c1post > 0.5) = 1;
            V_c1c2_pdf = bsxfun(@times,w1_unity, int_pdf) + bsxfun(@times,(1-w1_unity), V_seg_pdf);
            A_c1c2_pdf = bsxfun(@times,w1_unity, int_pdf) + bsxfun(@times,(1-w1_unity), A_seg_pdf);
        case 3 % probabilistic fusion, reweighting (covers possibility of always integrate/always segregate)
            fixed_weights = repmat(p_common,1,length(xrange),length(xrange)); %weights are fixed at the prior probability of one cause trials (fit to data) %JM might make this free parameter
            V_c1c2_pdf = bsxfun(@times,fixed_weights, int_pdf) + bsxfun(@times,(1-fixed_weights), V_seg_pdf);
            A_c1c2_pdf = bsxfun(@times,fixed_weights, int_pdf) + bsxfun(@times,(1-fixed_weights), A_seg_pdf);
        case 4 % probability matching
            
            
    end
end

%% Marginalize over internal variables using numerical integration
%integrate over xa xv values for every condition in condition vector
% try both using qtrapz adapted from acerbi

xpdf_V = bsxfun_normpdf(xrange_V, conds_V,V_sig);
xpdf_V = bsxfun(@rdivide, xpdf_V, qtrapz(xpdf_V, 2)); % Not multiplying by volume element
xpdf_A = bsxfun_normpdf(xrange_A, conds_A, A_sig);
xpdf_A = bsxfun(@rdivide, xpdf_A, qtrapz(xpdf_A, 3));  % Not multiplying by volume element(JM?)

if unity_judge
    prmat_unity = zeros(numel(conds_V), 2);
    %prmat_unity(:,1) = VestBMS_finalqtrapz(xpdf_V,xpdf_A,w1_unity);    % Not multiplying by volume element (xpdfs did not)
    prmat_unity(:,1) = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), w1_unity), 2), 3); %this is a recreation of what VestBMS_finalqtrapz does but not done in C.
    prmat_unity(:,2) = 1 - prmat_unity(:,1);
    % Fix probabilities
    prmat_unity = min(max(prmat_unity,0),1);
    %add in chance for random response, lambda
    prmat_unity = lambda/2 + (1-lambda)*prmat_unity;
end
if location_estimate %this is the location estimate when cause is unknown. Not sure if right.
    prmat_est_V = zeros(numel(conds_V), length(xrange));
    prmat_est_V(:,:) = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), V_c1c2_pdf), 2), 3); %this is a recreation of what VestBMS_finalqtrapz does but not done in C.
    prmat_est_A = zeros(numel(conds_A), length(xrange));
    prmat_est_A(:,:) = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), A_c1c2_pdf), 2), 3); %this is a recreation of what VestBMS_finalqtrapz does but not done in C.
    % Fix probabilities
    prmat_est_V= min(max(prmat_est_V,0),1);
    prmat_est_A= min(max(prmat_est_A,0),1);
    
    %add in chance for random response, lambda
    prmat_est_V = lambda/length(xrange) + (1-lambda)*prmat_est_V;
    prmat_est_A = lambda/length(xrange) + (1-lambda)*prmat_est_A;
    %lamda is the join probability
    prmat_est_A_reshape(:,1,1:size(prmat_est_A,2)) = prmat_est_A;
    prmat_sac = bsxfun(@times,prmat_est_V,prmat_est_A_reshape);
end


%% calculate negative log likelihood of data
if unity_judge && location_estimate %when fitting jointly
    nll_u = -1*sum(responses{1}.*log(prmat_unity));
    nll_u = sum(nll_u); %sum across bins
    nll_l = -1*sum(responses{2}.*log(prmat_sac));
    nll_l = sum(nll_l(:)); 
    
    nll = nll_u + nll_l; %total nll is sum of independent nll
    prmat{1} = prmat_unity;
    prmat{2} = prmat_sac;
elseif unity_judge
    nll = -1*sum(responses.*log(prmat_unity));
    nll = sum(nll(:)); %sum across bins
    prmat = prmat_unity;
    
elseif location_estimate %this is for the location estimate when cause is unknown
    nll = -1*sum(responses.*log(prmat_sac));
    nll = sum(nll(:)); %sum across bins
    prmat = prmat_sac;
end

end
