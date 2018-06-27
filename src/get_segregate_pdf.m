%% return pdf for full segregation condition
%
% -------------------
% Jeff Mohl
% 6/22/18
% -------------------
%
% Description: returns pdf assuming full segregation of A and V stimuli.
%
% Inputs:
% A_mu, V_mu, prior_mu: means of A, V, and prior components
% A_sig, V_sig, prior_sig: std  of same
% eval_pts: values for which the pdf should be calculated. Should either be
% a range (i.e. -100:100) or list of saccade endpoints if calculating
% likelihood

function [seg_pdf, A_pred,V_pred] = get_segregate_pdf(A_mu,V_mu,prior_mu,A_sig,V_sig,prior_sig,eval_pts)

        As_mu = (A_mu/A_sig^2 + prior_mu/prior_sig^2) / (1/A_sig^2 + 1/prior_sig^2); %incorporating prior sig on unisensory estimates
        Vs_mu = (V_mu/V_sig^2 + prior_mu/prior_sig^2) / (1/V_sig^2 + 1/prior_sig^2);
        adj_A_sig = sqrt((1/A_sig^2 + 1/prior_sig^2)^-1);
        adj_V_sig = sqrt((1/V_sig^2 + 1/prior_sig^2)^-1);
        A_pred = normpdf(eval_pts,As_mu,adj_A_sig); %evaluate likelihood at each actual saccade location
        V_pred = normpdf(eval_pts,Vs_mu,adj_V_sig);
        seg_pdf = (A_pred+V_pred)/2; %normalizing so that total probability sums to 1
end