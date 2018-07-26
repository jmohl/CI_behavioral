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

function n_loglikelihood = get_nll_CI_fixedAV(tar_pairs,endpoint_data,fixed_params,free_params)

%% extracting parameters for more readable code
prior_mu = fixed_params.prior_mu;

prior_sig = free_params(1);%sigma of centrality prior
prior_common = free_params(2);%prior on common cause

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
    
    %get fixed values from unimodal A and V for given target pair
    V_sig = fixed_params.V_sig(fixed_params.V_tars == V_tar); 
    A_sig = fixed_params.A_sig(fixed_params.A_tars == A_tar); 
    V_mu = fixed_params.V_mu(fixed_params.V_tars == V_tar); 
    A_mu = fixed_params.A_mu(fixed_params.A_tars == A_tar); 

    %xa and xv are assumed to be the mean of the unimodal responses, rather
    %than the true target locations
    xa = A_mu;
    xv = V_mu;
    
    %get pdf assuming integration
    int_pdf = get_integrate_pdf(xa,xv,prior_mu,A_sig,V_sig,prior_sig,endpoints);
    
    %get pdf assuming segregation
    seg_pdf = get_segregate_pdf(xa,xv,prior_mu,A_sig,V_sig,prior_sig,endpoints);
    
    %get posterior on common cause
    post_common = get_post_common(xa,xv,prior_mu,A_sig,V_sig,prior_sig,prior_common);
    
    %combine pdfs, weighted by post_common
    %note this is one 'optimal' way to do the task
    CI_pdf = post_common*int_pdf + (1-post_common)*seg_pdf;
    
    %calculate likelihood for observed saccades in this data by summing up
    %the values of the pdf and taking the negative log of that normalized
    %by the total number of saccades. This is the part I'm least sure is
    %correct.
    nll_array(i) = -1*sum(log(CI_pdf));
end

n_loglikelihood = sum(nll_array); %summing over all condition types.