%% get negative log likelihood for CI model
%
% -------------------
% Jeff Mohl
% 9/12/18
% -------------------
%
% Description: this function will calculate the negative log likelihood
% under the causal inference model.This code is adapted from get_nll_CI4 in
% the workspace code.
%
% Updated 9/12/18 to work using tidy_data structure instead of cell arrays,
% made it much much slower.

function n_loglikelihood = get_nll_CI(AV_pairs,fixed_params,free_params)
profile on
%% extracting parameters for more readable code
prior_mu = fixed_params.prior_mu;

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
nll_array = zeros(height(data),1);
for this_pair = 1:length(AV_pairs)
    A_tar = AV_pairs(this_pair,1);
    V_tar = AV_pairs(this_pair,2);
    this_data = data(data.A_tar == A_tar & data.V_tar == V_tar,:);
    endpoints = vertcat(this_data.valid_endpoints{:,1});
    endpoints = endpoints(:,1); %only taking x coord of endpoints, first column
    
    if (abs(A_tar) > 12)
        A_sig = Af_sig;         %switch which auditory sigma is used, based on target location
    else
        A_sig = Ac_sig;
    end
    
    % for this model, xa and xv are assumed to be centered on the true
    % target location. An option would be to have these be centered on the
    % unimodal distributions but that has resulted in some model
    % instability
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
    nll_array(this_pair) = -1*sum(log(CI_pdf));
end

n_loglikelihood = sum(nll_array); %summing over all condition types.
profile viewer