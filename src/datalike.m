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
global fitoptions %not sure if this slows things down or not
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

prior_type = model(4);
%give fit parameters sensible names
%minimum parameter set
V_sig = theta(1);%vis target sigma
A_sig = theta(2);%close aud target sigma
prior_sig = theta(3);%sigma of centrality prior
p_common = theta(4);%prior on common cause
lambda_uni = theta(5); %lapse probability
lambda_loc = theta(6);
prior_mu = 0; %fixed for now

%optional parameters
if ~prior_type == 1 %if non-normal prior is used
    resp_mu_A = theta(7);
    resp_mu_V = theta(8);
    resp_sig_A = theta(9);
    resp_sig_V = theta(10);
end

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
    sA(:,1) = xrange;
    sV(:,1) = xrange;
else
    %for multisensory trials response distributions are 2 dimensional
    %(joint A and V saccades)
    sA(:,1) = xrange;
    sV(:,1) = xrange;
    xrange_V(1,:,1) = xrange;
    xrange_A(1,1,:) = xrange;
end

%% (1) get sensory likelihood for visual and auditory p(xa|sa)
if prior_type ~= 1
    %normal prior makes this step irrelevant 
    likelihood_sA = bsxfun_normpdf(xrange_A,sA,A_sig);
    likelihood_sV = bsxfun_normpdf(xrange_V,sV,V_sig);
end
%% (2) get sensory priors 
if prior_type ~= 1
    if unisensory_loc
        locations = fitoptions.priorloc;
        prior_sA = get_unisensory_prior(prior_type,xrange,locations,prior_sig);
        prior_sV = get_unisensory_prior(prior_type,xrange,locations,prior_sig);
    else
        prior_c1 = get_joint_prior();
        prior_c2 = get_joint_prior();
    end
    
end

%% (3) find the posterior distribution for C = 1 case for all values of xa and xv;

switch CI_type 
    case 1 % bayesian causal inference posterior
        if prior_type == 1
            %c1 post is analytically solvable when prior is normal
            c1post = get_c1post(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,p_common);
        else
            %for all cases with non-normal priors, need to use numerical
            %methods
           
        end
end

%% (4) get posterior p(sA,sV|xa,xv,C)

if location_estimate
    if prior_type == 1
        %analytic solutions available for both of these cases, assuming a
        %normal prior
        int_pdf = get_integrate_pdf(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,xrange); %int pdf = 1x(xA)x(xV)x(eval_range) array. so for a given value of xA and xv, pdf is in 4th dim
        [~,A_seg_pdf,V_seg_pdf] = get_segregate_pdf(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,xrange); %seg pdf is 1x1x(xA)x(eval_range) array, so pdf is in 4th dimension for a given value of xA
    else
        
    end
    % combine likelihoods weighted by posterior on common cause
    %think this is more appropriately moved to the end of hte decision rule
    %step
end

if unisensory_loc
    if prior_type == 1
        [~,A_seg_pdf,V_seg_pdf] = get_segregate_pdf(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,xrange); %seg pdf is 1x1x(xA)x(eval_range) array, so pdf is in 4th dimension for a given value of xA
        post_sA = A_seg_pdf;
        post_sV = V_seg_pdf;
    else
       post_sA =likelihood_sA.*prior_sA';
       post_sA = post_sA./(1/sum(post_sA,'all')); 
       post_sV =likelihood_sV.*prior_sV';
       post_sV = post_sV./(1/sum(post_sV,'all')); 
    end  
end
%% (5) Incorporate decision rule
if location_estimate
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
if unity_judge
    switch CI_type
        case 1
            %judgement rule, if post > 0.5, choose unity
            w1_unity = zeros(size(c1post));
            w1_unity(c1post > 0.5) = 1;
            w1_unity(c1post == 0.5) = 0.5;
    case 2
        c1post = repmat(p_common,1,length(xrange),length(xrange));
        if unity_judge
            %if probabilistic fusion rule, judgement is fixed
            w1_unity = c1post;
        end
    end
end
    
%% (6) Marginalize over internal variables using numerical integration, to get estimates in terms of target locations
% integrate over xa xv values for every condition in condition vector
if ~unisensory_loc
    xpdf_V = bsxfun_normpdf(xrange_V, conds_V,V_sig);
    xpdf_V = bsxfun(@rdivide, xpdf_V, qtrapz(xpdf_V, 2)); % Not multiplying by volume element, this doesn't matter if the grid volumn = 1
    xpdf_A = bsxfun_normpdf(xrange_A, conds_A, A_sig);
    xpdf_A = bsxfun(@rdivide, xpdf_A, qtrapz(xpdf_A, 3));  % Not multiplying by volume element
else
    xpdf_V = bsxfun_normpdf(xrange_V, conds_V,V_sig);
    xpdf_V = bsxfun(@rdivide, xpdf_V, qtrapz(xpdf_V, 2));
    xpdf_A = bsxfun_normpdf(xrange_A, conds_A, A_sig);
    xpdf_A = bsxfun(@rdivide, xpdf_A, qtrapz(xpdf_A, 2)); 
end

%unisensory localization fits
if unisensory_loc 
    post_sA_rescale(1,:,:) = squeeze(post_sA);
    post_sV_rescale(1,:,:) = squeeze(post_sV);    

    prmat_est_A = zeros(numel(conds_A), length(xrange));
    prmat_est_A(:,:) = qtrapz(bsxfun(@times, xpdf_A, post_sA_rescale), 2);
    prmat_est_V = zeros(numel(conds_V), length(xrange));
    prmat_est_V(:,:) = qtrapz(bsxfun(@times, xpdf_V, post_sV_rescale), 2);
    %fix probabilities
    prmat_est_V= rdivide(min(max(prmat_est_V,0),1),sum(prmat_est_V,2));%division here rescales so sums to 1.
    prmat_est_A= rdivide(min(max(prmat_est_A,0),1),sum(prmat_est_A,2));
%     figure % for debugging
%     hold on
%     plot(xrange,prmat_est_A(1,:))
%     plot(xrange,prmat_est_A(2,:))
%     plot(xrange,prmat_est_A(3,:))
%     plot(xrange,prmat_est_A(4,:))
%     xticks(conds_A)
%     grid on
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
    if normal_prior %JM todo remove
        resp_prior_A = .5 * normpdf(xrange,-resp_mu_A,resp_sig_A) + .5* normpdf(xrange,resp_mu_A,resp_sig_A);
        resp_prior_V = .5 * normpdf(xrange,-resp_mu_V,resp_sig_A) + .5* normpdf(xrange,resp_mu_V,resp_sig_V);
        resp_A =A_c2_pdf .* resp_prior_A; % doing this makes the probabilities really small, so might need to rescale by the total so it sums to 1.
        resp_V =V_c2_pdf .* resp_prior_V;
        resp_AV =  AV_c1_pdf .* (0.5*(resp_prior_A) + 0.5*(resp_prior_V)); %TODO probably not this
    else
        resp_A =A_c2_pdf;
        resp_V =V_c2_pdf;
        resp_AV = AV_c1_pdf;
    end
    
    prmat_est_AV_c1 = zeros(numel(conds_V), length(xrange));
    prmat_est_AV_c1(:,:) = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), resp_AV), 2), 3); %Estimate for single saccade conditions, same for sA and sV, so labeling AV
    
    prmat_est_V_c2 = zeros(numel(conds_V), length(xrange));
    prmat_est_V_c2(:,:) = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), resp_V), 2), 3); %this is a recreation of what VestBMS_finalqtrapz does but not done in C.
    
    prmat_est_A_c2 = zeros(numel(conds_A), length(xrange));
    prmat_est_A_c2(:,:) = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), resp_A), 2), 3); %this is a recreation of what VestBMS_finalqtrapz does but not done in C.
    % Fix probabilities
    %JM todo make sure to rescale here, currently am not
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
%% (7) calculate negative log likelihood of data under this posterior
if unisensory_loc
    nll_A = -1 * sum(responses{1}.*log(prmat_est_A),'all');
    nll_V = -1 * sum(responses{2}.*log(prmat_est_V),'all');
    nll = nll_A + nll_V;
    prmat{1} = prmat_est_A;
    prmat{2} = prmat_est_V;

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
