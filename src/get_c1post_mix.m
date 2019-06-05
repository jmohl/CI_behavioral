%% calculate posterior probability of common cause for mixture prior
%
% -------------------
% Jeff Mohl - in progress
% 5/30/19
% -------------------
%
% Description: if inputs xa and xv are scalar, returns scalar
% probability of common cause between 0 and 1. If mu inputs are vectors,
% returns an array of probabilities for every combination of xa, xv. This
% is for the analytically solvable case where the prior is a single normal
% distribution.
%
% Inputs:
% Likelihood and prior distributions, prior probability of c=1 (p_common)
% Outputs:
% posterior distribution on c=1 for 1x(xv)x(xa) matrix

function c1post = get_c1post_mix(likelihood_xA, likelihood_xV,prior,p_common)
%c1post = p(xa,xv|C=1)*pr(C=1)/[p(xa,xv|C=1) * pr(C=1) + p(xa,xv|C=2)*(1-pr(C=1))]


%get p(xa,xv|C=1) * pr(C=1) [numerator of bayes rules
rs_prior(1,1,1,:) = prior; %need to put prior into the correct dimesion, which is the 4th dimension for the likelihoods
p_c1 = qtrapz(bsxfun(@times,bsxfun(@times,likelihood_xA,likelihood_xV),rs_prior),4)*p_common;

% get p(xa,xv|c=2) * (1-pr(c=1))
p_A_c2 = qtrapz(bsxfun(@times,likelihood_xA,rs_prior),4)*(1-p_common);
p_V_c2 = qtrapz(bsxfun(@times,likelihood_xV,rs_prior),4)*(1-p_common);
p_c2 = bsxfun(@times,p_A_c2,p_V_c2);

%final posterior combines the above computations
c1post = p_c1./(p_c1 + p_c2);

end
