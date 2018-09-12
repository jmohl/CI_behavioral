%% get negative log likelihood for CI model
%
% -------------------
% Jeff Mohl
% 5/30/18
% -------------------
%
% Description: this function will calculate the negative log likelihood
% under the causal inference model.This code is adapted from get_nll_CI4 in
% the workspace code.

function n_loglikelihood = get_nll_CI(endpoint_data,fixed_params,free_params)

%% extracting parameters for more readable code
prior_mu = fixed_params.prior_mu;
tar_pairs = fixed_params.AV_pairs;

V_sig = free_params(1);%.V_sig; %vis target sigma
Ac_sig = free_params(2);%.Ac_sig; %close aud target sigma
Af_sig = free_params(3);%.Af_sig; %far aud target sigma
prior_sig = free_params(4);%.prior_sig;%sigma of centrality prior
prior_common = free_params(5);%.prior_common;%prior on common cause

if prior_common > 1 || prior_common < 0 %can't have probability above 1 or below 0
    n_loglikelihood = 100000000;
    return
end

%% Calculate nll for each target pair and sum

nll_array = zeros(length(tar_pairs),1);
for i = 1:length(tar_pairs)
    A_tar = tar_pairs(i,1);
    V_tar = tar_pairs(i,2);
    endpoints = endpoint_data{i};
    
    if (abs(A_tar) > 12)
        A_sig = Af_sig;         %switch which auditory sigma is used, based on target location
    else
        A_sig = Ac_sig;
    end
    
    % currently xa and xv are assumed to be the true target location,
    % however this might change in the future if I decide to use the mean
    % of the reported location instead. this also might be incorrect, and
    % could be a sample taken from the mean instead? Could also make these
    % a free parameter to allow for bias.
    xa = A_tar;
    xv = V_tar;
    
    %get pdf assuming integration
    int_pdf = get_integrate_pdf(xa,xv,prior_mu,A_sig,V_sig,prior_sig,endpoints);
    
    %get pdf assuming segregation
    seg_pdf = get_segregate_pdf(xa,xv,prior_mu,A_sig,V_sig,prior_sig,endpoints);
    
    %get posterior on common cause
    post_common = get_post_common(xa,xv,prior_mu,A_sig,V_sig,prior_sig,prior_common);
    
    %combine pdfs, weighted by post_common
    %note this is one 'optimal' way to do the task
    CI_pdf = post_common*int_pdf + (1-post_common)*seg_pdf;
    
    %calculate likelihood for observed saccades in this data
    nll_array(i) = -1*sum(log(CI_pdf));
end

n_loglikelihood = sum(nll_array); %summing over all condition types.