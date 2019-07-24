%% Likelihood calculation for various models
%
% -------------------
% Jeff Mohl
% 1/30/19
% -------------------
%
% Description: this the updated code for calculating the log likelihood of
% a given dataset under each model.
%
% Inputs:
% conditions: cell array of conditions used with format {Atar1,Vtar1;
% Atar2,Vtar2;...}

%responses: cell array of responses, different for each task type
% unity: (nconditions)x2 array of responses for unity (column 1) vs
% separate (column 2)
% location: (nconditions)x101x101 array of sacade counts in 1 degree bins,
% where A saccade and V saccade locations define a point on the 101x101
% grid, with single saccade trials on the diagonal
% joint: cell array containing both of the above.

% Theta: parameters used for model
% Included in all models:
% theta(1) = aud sigma
% theta(2) = vis sigma
% theta(3) = prior sigma, variance of prior (not used in discrete prior case)
% theta(4) = p_common, prior probability of common cause
% theta(5) = lapse probability on unity judgement
% Optional thetas
% theta(6) = prior_mu2 (when using mixture of normals prior)
% theta(6) = prior_sig2 (when using mixture of normals prior)

% model: definition of model used
% model(1) = CI type: Bayes (1) or probabilistic fusion (2)
% model(2) = stimulus fusion: Bayes reweight(1), model selection(2), probabilistic fusion (3)
% model(3) = task/fit type: unity judgement (1), localization (2), joint fit (3), unisensory localization (4)
% model(4) = prior type: naive normal (1), discrete empirical (2), normal mixture empirical (3)

% eval_midpoint: evaluation points for numerical methods. Must match with
% saccade bins.

%% start of likelihood code
function [nll,prmat] = datalike(conditions,responses,theta,model,eval_midpoints)
%set up options 
debug = 0;
incorporate_unity_lapse = 1; %testing this feature, where the subject might randomly make 1 saccade on their localization trials for instance. Currently don't like.

CI_type = model(1);
combination_rule = model(2); %rule for combining sensory inputs, based on causal judgement 
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
A_sig = theta(1);%close aud target sigma
V_sig = theta(2);%vis target sigma
prior_sig = theta(3);%sigma of centrality prior
p_common = theta(4);%prior on common cause
lambda_uni = theta(5); %lapse probability
prior_mu = 0;

%optional parameters
if prior_type == 3 %if mixture prior used
    prior_mu2 = theta(6);
    prior_sig2 = theta(7);
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
if min(theta(1:3)) <= 0.1 || max(theta(4:5)) > 1 || min(theta(4:5)) < 0
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
    %for subsequent steps, format will always be
    %(conditions)x(xv)x(xa)x(pdf) for the pdf matrices.
    xrange_A(1,1,:) = xrange;
    xrange_V(1,:,1) = xrange;
    sA(1,1,1,:) = xrange; %ideally I would use different dimensions for this, but the matrices become too large
    sV(1,1,1,:) = xrange;

end

%test params
if debug
    A_sig = 5;
    V_sig = 2;
    prior_sig = 10;
    prior_mu = 15;
    p_common = 0.5;
end

%% (1) get sensory likelihood for visual and auditory p(xa|sa)
if prior_type ~= 1 %naive normal prior is handled analytically in subsequent steps
    likelihood_xA = bsxfun_normpdf(xrange_A,sA,A_sig);
    likelihood_xV = bsxfun_normpdf(xrange_V,sV,V_sig);
end
%% (2) get stimulus location priors 
if prior_type ~= 1 %naive normal prior is handled analytically in subsequent steps
    switch prior_type
        case 1
            locations = 0;
        case 2
            locations = [-24 -18 -12 -6 6 12 18 24]; %use actual target locations
        case 3
            locations = [prior_mu2,prior_mu,-prior_mu2];%symmetric around 0
            prior_sig = [prior_sig2,prior_sig,prior_sig2];
    end
    prior = get_unisensory_prior(prior_type,xrange,locations,prior_sig);
end

%% (3) find the posterior distribution for C = 1 case for all values of xa and xv;

if prior_type == 1
    %c1 post is analytically solvable when prior is normal
    c1post = get_c1post(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,p_common); %switched prior mu with 0 here (3rd param) because I'm being lazy, will remove once I add the new CI priors.
else
    %for all cases with non-normal priors, need to use integration
    c1post = get_c1post_mix(likelihood_xA,likelihood_xV, prior,p_common);
end

%% (4) get p(sA,sV|xa,xv,C) if estimating location

if location_estimate
    if prior_type == 1
        %analytic solutions available for both of these cases, assuming a
        %normal prior
        int_pdf = get_integrate_pdf(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,xrange); %int pdf = 1x(xA)x(xV)x(eval_range) array. so for a given value of xA and xv, pdf is in 4th dim
        [~,A_seg_pdf,V_seg_pdf] = get_segregate_pdf(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,xrange); %seg pdf is 1x1x(xA)x(eval_range) array, so pdf is in 4th dimension for a given value of xA
    else
        %c1 case
        %prior needs to go in the 4th dimensio
        rs_prior(1,1,1,:) = prior; 
        int_pdf = bsxfun(@times,bsxfun(@times,likelihood_xV,likelihood_xA),rs_prior);
        int_pdf = bsxfun(@rdivide,int_pdf,sum(int_pdf,4)); %can't do it this way because it ruins things at the edges. TODO fix this
        %c2 case
        A_seg_pdf = bsxfun(@times,likelihood_xA, rs_prior);
        A_seg_pdf = bsxfun(@rdivide,A_seg_pdf,bsxfun(@rdivide,sum(A_seg_pdf,4),sum(likelihood_xA,4)*sum(prior)));
        V_seg_pdf = bsxfun(@times,likelihood_xV, rs_prior);
        V_seg_pdf = bsxfun(@rdivide,V_seg_pdf,bsxfun(@rdivide,sum(V_seg_pdf,4),sum(likelihood_xV,4)*sum(prior)));
    end
end

if unisensory_loc
    if prior_type == 1
        [~,A_seg_pdf,V_seg_pdf] = get_segregate_pdf(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,xrange); %seg pdf is 1x1x(xA)x(eval_range) array, so pdf is in 4th dimension for a given value of xA
        post_sA = A_seg_pdf;
        post_sV = V_seg_pdf;
    else
        %prior needs to go in first dimension
       post_sA =likelihood_xA.*prior';
       post_sA = post_sA./(1/sum(post_sA,'all')); 
       post_sV =likelihood_xV.*prior';
       post_sV = post_sV./(1/sum(post_sV,'all')); 
    end  
end
%% (5) Incorporate decision rule
if location_estimate 
    switch combination_rule
        case 1 %if bayesian reweighting
            %reweight by posterior
            if incorporate_unity_lapse %jm todo remove
                w_unity = c1post*(1-lambda_uni) + lambda_uni/2; %there is the possibility of making a mistake about number of saccades, that will propogate to the location report
                w_unity_lapse = c1post*(1-lambda_uni);
            else
                w_unity = c1post;
            end
        case 2 %if model selection, choose the best model and use that one exclusively (set weight on alternative to 0)
            w_unity = zeros(size(c1post));
            if incorporate_unity_lapse %JM todo remove
                w_unity(c1post > 0.5) = 1-lambda_uni/2; %there is the possibility of making a mistake about number of saccades, that will propogate to the location report
                w_unity_lapse(c1post > 0.5) = 1-lambda_uni;
            else
                w_unity(c1post > 0.5) = 1;
            end
        case 3 % probabilistic fusion, reweighting (covers possibility of always integrate/always segregate)
            w_unity = repmat(p_common,1,length(xrange),length(xrange)); %weights are fixed at the prior probability of one cause trials (fit to data) %JM might make this free parameter
        case 4 % probability matching 0 - not actually sure if this is differentiable from reweighting with how I've done the task,
            %todo
    end

    %if option enabled, compensate for the possibility of randomly making 1
    %saccade by including a scaled projection
    if incorporate_unity_lapse
        %c2 pdfs are handled the same, simply scaling down by the
        %probability of randomly making only a single saccade
        V_c2_pdf = bsxfun(@times,(1-w_unity), V_seg_pdf);
        A_c2_pdf = bsxfun(@times,(1-w_unity), A_seg_pdf);
        %c1 probability is more complicated, because the lapse saccade
        %could either be taken from the c1, A, or V distributions. Here I
        %am assuming that the lapses are coming equally from the A and V
        %distributions.
        AV_c1_pdf = bsxfun(@times, w_unity_lapse, int_pdf);
        lapse_c1_pdf = (A_seg_pdf + V_seg_pdf)*lambda_uni/4; 
        AV_c1_pdf = AV_c1_pdf + lapse_c1_pdf;
        resp_A_c2 = A_c2_pdf;
        resp_V_c2 = V_c2_pdf;
        resp_AV_c1 = AV_c1_pdf;
    else
        AV_c1_pdf = bsxfun(@times,w_unity, int_pdf); %splitting these because my responses are split between 1 and 2 saccade cases.
        V_c2_pdf = bsxfun(@times,(1-w_unity),V_seg_pdf);
        A_c2_pdf = bsxfun(@times,(1-w_unity), A_seg_pdf);
        resp_A_c2 = A_c2_pdf;
        resp_V_c2 = V_c2_pdf;
        resp_AV_c1 = AV_c1_pdf;
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
    %decision rule includes chance for random choice
    resp_unity = lambda_uni/2 + (1-lambda_uni)*w1_unity;
end

if unisensory_loc
    %assume chance for random saccades, but that's it
    resp_sA = post_sA;
    resp_sV = post_sV;
end
%% (6) Marginalize over internal variables using numerical integration, to get estimates in terms of target locations
% integrate over xa xv values for every condition in condition vector
if ~unisensory_loc
    xpdf_V = bsxfun_normpdf(xrange_V, conds_V,V_sig);
    xpdf_V = bsxfun(@rdivide, xpdf_V, qtrapz(xpdf_V, 2)); % Not multiplying by volume element, this doesn't matter if the grid volumn = 1
    xpdf_A = bsxfun_normpdf(xrange_A, conds_A, A_sig);
    xpdf_A = bsxfun(@rdivide, xpdf_A, qtrapz(xpdf_A, 3));  % Not multiplying by volume element
else %unisensory loc just uses 3d grids
    xpdf_V = bsxfun_normpdf(xrange_V, conds_V,V_sig);
    xpdf_V = bsxfun(@rdivide, xpdf_V, qtrapz(xpdf_V, 2));
    xpdf_A = bsxfun_normpdf(xrange_A, conds_A, A_sig);
    xpdf_A = bsxfun(@rdivide, xpdf_A, qtrapz(xpdf_A, 2)); 
end

%unisensory localization fits
if unisensory_loc 
    resp_sA_rescale(1,:,:) = squeeze(resp_sA);
    resp_sV_rescale(1,:,:) = squeeze(resp_sV);    

    prmat_est_A = zeros(numel(conds_A), length(xrange));
    prmat_est_A(:,:) = qtrapz(bsxfun(@times, xpdf_A, resp_sA_rescale), 2);
    prmat_est_V = zeros(numel(conds_V), length(xrange));
    prmat_est_V(:,:) = qtrapz(bsxfun(@times, xpdf_V, resp_sV_rescale), 2);
    %fix probabilities
    prmat_est_V= rdivide(min(max(prmat_est_V,0),1),sum(prmat_est_V,2));%division here rescales so sums to 1.
    prmat_est_A= rdivide(min(max(prmat_est_A,0),1),sum(prmat_est_A,2));
end
    
%estimate percent of single saccade trials, integrate per condition
if unity_judge
    prmat_unity = zeros(numel(conds_V), 2);
    prmat_unity(:,1) = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), resp_unity), 2), 3); %this is a recreation of what VestBMS_finalqtrapz does but not done in C.
    prmat_unity(:,2) = 1 - prmat_unity(:,1);
    % Fix probabilities
    prmat_unity = min(max(prmat_unity,0),1);
end

%estimate saccade endpoint locations, integrate per condition
if location_estimate %this is the location estimate when cause is unknown. Not sure if right.
    prmat_AV_c1 = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), resp_AV_c1), 2), 3); %this is a recreation of what VestBMS_finalqtrapz does but not done in C.
    prmat_A_c2 = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), resp_A_c2), 2), 3);
    prmat_V_c2 = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), resp_V_c2), 2), 3);
    prmat_V_shape(:,:,1) = prmat_V_c2;
    prmat_A_shape(:,1,:) = prmat_A_c2;
    prmat_AV_c2 = bsxfun(@times,prmat_A_shape,prmat_V_shape);%problem here is normalizing to make the probabilities match, because this is going to significantly reduce them
    %renormalize so that total probability will sum to 1 for each condition
    %after adding c=1 case
    prmat_AV_c2 = bsxfun(@rdivide,prmat_AV_c2,sum(prmat_V_c2,4));
    
    %diagonalize 
    diag_array = diag(ones(length(xrange),1));
    diag_array = logical(diag_array);
    prmat_sac = prmat_AV_c2;
    %add C=1 case to diagonal, producing final cxSvxSa matrix
    prmat_sac(:,diag_array) = prmat_sac(:,diag_array) + squeeze(prmat_AV_c1(:,1,1,:));
    
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
