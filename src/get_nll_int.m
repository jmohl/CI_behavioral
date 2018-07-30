%% get negative log likelihood for separated model
%
% -------------------
% Jeff Mohl
% 7/30/18
% -------------------
%
% Description: this function will calculate the negative log likelihood
% under a model which always assumes separate targets. 

function n_loglikelihood = get_nll_int(tar_pairs,endpoint_data,fixed_params,free_params)

%% extracting parameters for more readable code
prior_mu = fixed_params.prior_mu;

V_sig = free_params(1);%.V_sig; %vis target sigma
Ac_sig = free_params(2);%.Ac_sig; %close aud target sigma
Af_sig = free_params(3);%.Af_sig; %far aud target sigma
prior_sig = free_params(4);%.prior_sig;%sigma of centrality prior


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
    xa = A_tar;
    xv = V_tar;
    
    %get pdf assuming segregation
    %get pdf assuming integration
    int_pdf = get_integrate_pdf(xa,xv,prior_mu,A_sig,V_sig,prior_sig,endpoints);
    
    %calculate likelihood for observed saccades in this data
    nll_array(i) = -1*sum(log(int_pdf));
end

n_loglikelihood = sum(nll_array); %summing over all condition types.