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
% theta(5) lambda_uni: lapse rate, probability randomly making unity judgement
% theta(6) lambda_loc: lapse rate, probability of making a random saccade
% model(1) Causal inference type (1=bayesian,2=probabilistic fusion, 3 fixed criterion (todo))
% model(2) combination rule for localization(1=posterior reweighting, 2 = model selection, 3 = probabilistic fusion (todo, make fusion rate free param?) 4= probability matching (todo))
% model(3) Response type (1 = unity judgement, 2 = location, 3 = joint)
%


%% start of likelihood code - currently only working for unity judgement
function [nll,prmat] = datalike(conditions,responses,theta,model,eval_midpoints)
% Causal inference type
CI_type = model(1);
%localization strategy
combination_rule = model(2); %rule for combining sensory inputs, based on causal judgement 
%task type

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

%give fit parameters sensible names
V_sig = theta(1);%vis target sigma
A_sig = theta(2);%close aud target sigma
prior_sig = theta(3);%sigma of centrality prior
p_common = theta(4);%prior on common cause
lambda_uni = theta(5); %lapse probability
lambda_loc = theta(6);
prior_mu = 0; %fixed for now

if unisensory_loc
    conds_A = conditions{1};
    conds_V = conditions{2};
else
    conds_A = conditions(:,1);
    conds_V = conditions(:,2);
end

%HACK fminsearch does not allow bounds, so I'm including this really hacky
%way to require that the likelihood is not calculated for impossible
%values. In the future I might switch to bads or fminsearchbnd, both of
%which are third party but allow bounds. 
if min(theta(1:3)) <= 0.1 || max(theta(4:6)) > 1 || min(theta(4:6)) < 0
    nll = 1e10;
    return;
end
xrange = eval_midpoints; %range for integration.
if unisensory_loc
    xrange_V(1,:) = xrange;
    xrange_A(1,:) = xrange;
else
    %for multisensory trials response distributions are 2 dimensional
    %(joint A and V saccades)
    xrange_V(1,:,1) = xrange;
    xrange_A(1,1,:) = xrange;
end

%% find the posterior distribution for C = 1 case for all values of xa and xv;

switch CI_type 
    case 1 % bayesian causal inference posterior
        c1post = get_c1post(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,p_common);
        if unity_judge
               %judgement rule, if post > 0.5, choose unity
            w1_unity = zeros(size(c1post));
            w1_unity(c1post > 0.5) = 1;
            w1_unity(c1post == 0.5) = 0.5;
        end
    case 2
        c1post = repmat(p_common,1,length(xrange),length(xrange));
        if unity_judge
            %if probabilistic fusion rule, judgement is fixed
            w1_unity = c1post;
        end
end

%% get likelihood functions for A location and V location
if location_estimate || unisensory_loc
    int_pdf = get_integrate_pdf(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,xrange); %int pdf = 1x(xA)x(xV)x(eval_range) array. so for a given value of xA and xv, pdf is in 4th dim
    
    [~,A_seg_pdf,V_seg_pdf] = get_segregate_pdf(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,xrange); %seg pdf is 1x1x(xA)x(eval_range) array, so pdf is in 4th dimension for a given value of xA
    
    % combine likelihoods weighted by posterior on common cause
    switch combination_rule
        case 1 %if bayesian reweighting
            AV_c1_pdf = bsxfun(@times,c1post, int_pdf) ; %splitting these because my responses are split between 1 and 2 saccade cases.
            V_c2_pdf = bsxfun(@times,(1-c1post), V_seg_pdf);
            A_c2_pdf = bsxfun(@times,(1-c1post), A_seg_pdf);
        case 2 %if model selection, choose the best model and use that one exclusively (set weight on alternative to 0)
            w1_unity = zeros(size(c1post));
            w1_unity(c1post > 0.5) = 1;
            AV_c1_pdf = bsxfun(@times,w1_unity, int_pdf);
            V_c2_pdf = bsxfun(@times,(1-w1_unity), V_seg_pdf);
            A_c2_pdf = bsxfun(@times,(1-w1_unity), A_seg_pdf);
        case 3 % probabilistic fusion, reweighting (covers possibility of always integrate/always segregate)
            fixed_weights = repmat(p_common,1,length(xrange),length(xrange)); %weights are fixed at the prior probability of one cause trials (fit to data) %JM might make this free parameter
            AV_c1_pdf = bsxfun(@times,fixed_weights, int_pdf);
            V_c2_pdf = bsxfun(@times,(1-fixed_weights), V_seg_pdf);
            A_c2_pdf = bsxfun(@times,(1-fixed_weights), A_seg_pdf);
        case 4 % probability matching
            %todo
            
    end
end

%% Marginalize over internal variables using numerical integration
%integrate over xa xv values for every condition in condition vector
% try both using qtrapz adapted from acerbi
if ~unisensory_loc
    xpdf_V = bsxfun_normpdf(xrange_V, conds_V,V_sig);
    xpdf_V = bsxfun(@rdivide, xpdf_V, qtrapz(xpdf_V, 2)); % Not multiplying by volume element
    xpdf_A = bsxfun_normpdf(xrange_A, conds_A, A_sig);
    xpdf_A = bsxfun(@rdivide, xpdf_A, qtrapz(xpdf_A, 3));  % Not multiplying by volume element(JM?)
end
%estimate percent of single saccade trials, integrate per condition
if unity_judge
    prmat_unity = zeros(numel(conds_V), 2);
    %prmat_unity(:,1) = VestBMS_finalqtrapz(xpdf_V,xpdf_A,w1_unity);    % Not multiplying by volume element (xpdfs did not)
    prmat_unity(:,1) = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), w1_unity), 2), 3); %this is a recreation of what VestBMS_finalqtrapz does but not done in C.
    prmat_unity(:,2) = 1 - prmat_unity(:,1);
    % Fix probabilities
    prmat_unity = min(max(prmat_unity,0),1);
    %add in chance for random response, lambda
    prmat_unity = lambda_uni/2 + (1-lambda_uni)*prmat_unity;
end

%estimate saccade endpoint locations, integrate per condition
if location_estimate %this is the location estimate when cause is unknown. Not sure if right.
    prmat_est_AV_c1 = zeros(numel(conds_V), length(xrange));
    prmat_est_AV_c1(:,:) = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), AV_c1_pdf), 2), 3); %Estimate for single saccade conditions, same for sA and sV, so labeling AV
    
    prmat_est_V_c2 = zeros(numel(conds_V), length(xrange));
    prmat_est_V_c2(:,:) = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), V_c2_pdf), 2), 3); %this is a recreation of what VestBMS_finalqtrapz does but not done in C.
    
    prmat_est_A_c2 = zeros(numel(conds_A), length(xrange));
    prmat_est_A_c2(:,:) = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), A_c2_pdf), 2), 3); %this is a recreation of what VestBMS_finalqtrapz does but not done in C.
    % Fix probabilities
    prmat_est_V= min(max(prmat_est_V_c2,0),1);
    prmat_est_A= min(max(prmat_est_A_c2,0),1);
    prmat_est_AV = min(max(prmat_est_AV_c1,0),1);
    
    %add in chance for random response, lambda
    prmat_est_V = lambda_loc/length(xrange) + (1-lambda_loc)*prmat_est_V;
    prmat_est_A = lambda_loc/length(xrange) + (1-lambda_loc)*prmat_est_A;
    prmat_est_AV = lambda_loc/length(xrange) + (1-lambda_loc)*prmat_est_AV;

    prmat_est_A_reshape(:,1,1:size(prmat_est_A,2)) = prmat_est_A;
    prmat_sac_c2 = bsxfun(@times,prmat_est_V,prmat_est_A_reshape);
    %diagonalize prmat_est_AV, so can be combined with the c2 estimates
    %which are 2d
    prmat_est_AV_diag = zeros(size(prmat_sac_c2));
    for c=1:length(conditions)
        prmat_est_AV_diag(c,:,:) = diag(prmat_est_AV(c,:));
    end
    prmat_sac = prmat_sac_c2 + prmat_est_AV_diag;
end

if unisensory_loc %unisensory localization fits
    A_seg_pdf_rescale(1,:,:) = squeeze(A_seg_pdf);
    V_seg_pdf_rescale(1,:,:) = squeeze(V_seg_pdf);
    xpdf_V = bsxfun_normpdf(xrange_V, conds_V,V_sig);
    xpdf_V = bsxfun(@rdivide, xpdf_V, qtrapz(xpdf_V, 2));
    xpdf_A = bsxfun_normpdf(xrange_A, conds_A,A_sig);
    xpdf_A = bsxfun(@rdivide, xpdf_A, qtrapz(xpdf_A, 2));
    
    prmat_est_V_uni = zeros(numel(conds_V), length(xrange));
    prmat_est_V_uni(:,:) = qtrapz(bsxfun(@times, xpdf_V, V_seg_pdf_rescale), 2);
    prmat_est_A_uni = zeros(numel(conds_A), length(xrange));
    prmat_est_A_uni(:,:) = qtrapz(bsxfun(@times, xpdf_A, A_seg_pdf_rescale), 2);
    
    %fix probabilities
    prmat_est_V_uni= min(max(prmat_est_V_uni,0),1);
    prmat_est_A_uni= min(max(prmat_est_A_uni,0),1);
    
    %add lapse rate
    prmat_est_V_uni = lambda_loc/length(xrange) + (1-lambda_loc)*prmat_est_V_uni; %chance to make a lapse sacade is uniform across space. which is not technically true but an approximation
    prmat_est_A_uni = lambda_loc/length(xrange) + (1-lambda_loc)*prmat_est_A_uni;
    
end
    

%% calculate negative log likelihood of data
if unisensory_loc
    nll_A = -1 * sum(responses{1}.*log(prmat_est_A_uni),'all');
    nll_V = -1 * sum(responses{2}.*log(prmat_est_V_uni),'all');
    nll = nll_A + nll_V;
    prmat{1} = prmat_est_A_uni;
    prmat{2} = prmat_est_V_uni;

elseif unity_judge && location_estimate %when fitting jointly
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
