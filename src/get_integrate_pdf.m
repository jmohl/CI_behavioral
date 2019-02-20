%% return pdf for full integration condition
%
% -------------------
% Jeff Mohl
% 2/14/19
% -------------------
%
% Description: returns pdf assuming full integration of A and V stimuli. If
% xa and xv are scalar values a single pdf is returns for the range of
% values in eval_range. If xa and xv are vectors, will return a 1xMxNxD array
% where M = N = length(xa), and D = length(eval_range).
%
% Inputs:
% A_mu, V_mu, prior_mu: means of A, V, and prior components
% A_sig, V_sig, prior_sig: std  of same
% eval_pts: values for which the pdf should be calculated. Should either be
% a range (i.e. -100:100) or list of saccade endpoints if calculating
% likelihood

function [int_pdf, int_mu, int_sig] = get_integrate_pdf(xa,xv,prior_mu,A_sig,V_sig,prior_sig,eval_range)
if length(xa) == 1
    int_mu = (xv/V_sig^2 + xa/A_sig^2 + prior_mu/prior_sig^2)/(1/A_sig^2 + 1/V_sig^2 + 1/prior_sig^2);
    int_sig =sqrt((1/A_sig^2 + 1/V_sig^2 + 1/prior_sig^2)^-1);
    int_pdf = normpdf(eval_range,int_mu,int_sig);
else
    int_mu = (bsxfun(@plus,xv./V_sig^2, xa./A_sig^2) + prior_mu/prior_sig^2)/(1/A_sig^2 + 1/V_sig^2 + 1/prior_sig^2);
    int_sig =sqrt((1/A_sig^2 + 1/V_sig^2 + 1/prior_sig^2)^-1);
    expanded_range(1,1,1,:) = eval_range; %putting eval range in 4th dim to allow for integration later on.
    int_pdf = bsxfun_normpdf(expanded_range,int_mu,int_sig);
end

end