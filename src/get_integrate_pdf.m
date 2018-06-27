%% return pdf for full integration condition
%
% -------------------
% Jeff Mohl
% 6/22/18
% -------------------
%
% Description: returns pdf assuming full integration of A and V stimuli
%
% Inputs:
% A_mu, V_mu, prior_mu: means of A, V, and prior components
% A_sig, V_sig, prior_sig: std  of same
% eval_pts: values for which the pdf should be calculated. Should either be
% a range (i.e. -100:100) or list of saccade endpoints if calculating
% likelihood

function int_pdf = get_integrate_pdf(A_mu,V_mu,prior_mu,A_sig,V_sig,prior_sig,eval_pts)
        mu_pred = (V_mu/V_sig^2 + A_mu/A_sig^2 + prior_mu/prior_sig^2)/(1/A_sig^2 + 1/V_sig^2 + 1/prior_sig^2);
        sig_pred =sqrt((1/A_sig^2 + 1/V_sig^2 + 1/prior_sig^2)^-1);
        int_pdf = normpdf(eval_pts,mu_pred,sig_pred);
end