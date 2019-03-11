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
% model(1) Model type (1=Bayesian reweighting v1, 2=bayesian reweighting v2, 3=probabilistic fusion, 4 = model selection)
% model(2) Response type (1=unity judgement, 2 = location, 3 = fit joint)
% model(3) estimation proceedure (1=numerical integration, 2 = analytic)

%TODO: this is getting a little confusing with the factorial model types.
%need to double check that the proper values are being used at all the
%different stages for if statements.


%% start of likelihood code - currently only working for unity judgement
function [nll,prmat] = datalike(conditions,responses,theta,model,eval_midpoints)

CI_type = model(1);
unity_judge = model(2) == 1;
location_estimate = model(2) == 2;
if model(2) == 3 %do joint fit
    unity_judge = 1;
    location_estimate = 1;
end

method = model(3); %fit using numerical integration, or analytically


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

if CI_type ~= 3 % probabilistic fusion model, c1post is fixed at p_common
    c1post = get_c1post(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,p_common);
end
if unity_judge
    %judgement rule, if post > 0.5, choose unity
    if CI_type == 3
        w1_unity = repmat(p_common,1,length(xrange),length(xrange));
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
    if CI_type == 2 %if bayesian reweighting
        V_c1c2_pdf = bsxfun(@times,c1post, int_pdf) + bsxfun(@times,(1-c1post), V_seg_pdf);
        A_c1c2_pdf = bsxfun(@times,c1post, int_pdf) + bsxfun(@times,(1-c1post), A_seg_pdf);
    elseif CI_type == 4 %if model selection, choose the best model and use that one exclusively (set weight on alternative to 0)
        w1_unity = zeros(size(c1post));
        w1_unity(c1post > 0.5) = 1;
        V_c1c2_pdf = bsxfun(@times,w1_unity, int_pdf) + bsxfun(@times,(1-w1_unity), V_seg_pdf);
        A_c1c2_pdf = bsxfun(@times,w1_unity, int_pdf) + bsxfun(@times,(1-w1_unity), A_seg_pdf);
    end
end

%% Marginalize over internal variables using numerical integration
%integrate over xa xv values for every condition in condition vector
% try both using qtrapz adapted from acerbi

if method == 1 %method is numerical integration
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
        if CI_type == 1 %old bayesian reweighting method, saccades are unlabeled and all saccades are lumped together
            prmat_sac = (prmat_est_V + prmat_est_A)/2; %normalized
        else %new bayesian reweighting, A and V saccades are labeled, using joint distribution instead of sum of distributions
            prmat_est_A_reshape(:,1,1:size(prmat_est_A,2)) = prmat_est_A;
            prmat_sac = bsxfun(@times,prmat_est_V,prmat_est_A_reshape);
        end
        
    end
elseif method == 2 %analytic solution exists 3/6/19 don't think this actually works, should probably remove
    if unity_judge && CI_type == 3;
        prmat_unity = zeros(numel(conds_V), 2);
        prmat_unity(:,1) = repmat(p_common,size(prmat_unity(:,1)));
        prmat_unity(:,2) = 1 - prmat_unity(:,1);
        % Fix probabilities
        prmat_unity = min(max(prmat_unity,0),1);
        %add in chance for random response, lambda % changed my mind, for
        %this model the only parameter is p_common, which is basically
        %lambda
        %prmat_unity = lambda/2 + (1-lambda)*prmat_unity;
    end
    if location_estimate %analytic solution can be found assuming the integrate and segregate pdfs are independent and combined according to constant weight
        %for each condition, can solve for mu and sigma analytically. use
        %these to make the pr_mat
        [~, intx_mu,intx_sig] = get_integrate_pdf(conditions(:,1),conditions(:,2),prior_mu,A_sig,V_sig,prior_sig,xrange);
        %get params of marginal distribution
        int_mu = (bsxfun(@plus,bsxfun(@plus,conditions(:,2)./V_sig^2, conditions(:,1)./A_sig^2), intx_mu./intx_sig^2))/(1/A_sig^2 + 1/V_sig^2 + 1/intx_sig^2);
        int_sig =sqrt((1/A_sig^2 + 1/V_sig^2 + 1/intx_sig^2)^-1);
        int_pdf = bsxfun_normpdf(xrange,int_mu,int_sig);
        
        [~, ~,~,Asx_mu, Vsx_mu,Asx_sig,Vsx_sig] = get_segregate_pdf(conditions(:,1),conditions(:,2),prior_mu,A_sig,V_sig,prior_sig,xrange);
        As_mu = bsxfun(@plus,conditions(:,1)./A_sig^2, Asx_mu./Asx_sig^2)/(1/A_sig^2 +1/Asx_sig^2);
        As_sig = sqrt((1/A_sig^2 + 1/Asx_sig^2)^-1);
        Vs_mu = bsxfun(@plus,conditions(:,1)./A_sig^2, Vsx_mu./Vsx_sig^2)/(1/A_sig^2 +1/Vsx_sig^2);
        Vs_sig = sqrt((1/A_sig^2 + 1/Vsx_sig^2)^-1);
        seg_pdf = bsxfun(@times,bsxfun_normpdf(xrange_A,As_mu,As_sig),bsxfun_normpdf(xrange_V,Vs_mu,Vs_sig)) ;
        
        %combine = int_pdf added to diagonal of seg pdf, with both weighted
        %by post_common
        I_mat = [];
        I_mat(1,:,:) = eye(length(xrange),length(xrange));
        I_mat = repmat(I_mat,length(conditions),1,1);
        int_pdf_diag = bsxfun(@times,I_mat,int_pdf);
        
        prmat_sac =  bsxfun(@times,int_pdf_diag, c1post) + bsxfun(@times,seg_pdf,(1-c1post));
        %add random saccades
        prmat_sac = lambda/length(xrange)^2 + (1-lambda) * prmat_sac;
    end
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
