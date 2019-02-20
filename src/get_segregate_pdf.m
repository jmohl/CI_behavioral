%% return pdf for full segregation condition
%
% -------------------
% Jeff Mohl
% 2/14/19
% -------------------
%
% Description: returns pdf assuming full segregation of A and V stimuli.
% Returns the individual pdfs for the A and V targets, as well as a
% seg_pdf which assumes that both the A and V distributions contribute to
% the total saccade distribution equally (sums A and V distributions). If
% the inputs xa and xv are scalar, returns a single distribution along
% eval_range. If the xa and xv inputs are vectors, will return a 1xMxNxD array
% where M = N = length(xa), and D = length(eval_range).
%
% Inputs:
% A_mu, V_mu, prior_mu: means of A, V, and prior components
% A_sig, V_sig, prior_sig: std  of same
% eval_pts: values for which the pdf should be calculated. Should either be
% a range (i.e. -100:100) or list of saccade endpoints if calculating
% likelihood
%
% Outputs: 
% seg_pdf, the sum of A_pdf and V_pdf. A
% A_pdf, normal pdf of auditory target sa given xa and prior
% V_pdf, same but for visual percept xv.
function [seg_pdf, A_pdf,V_pdf, As_mu, Vs_mu,As_sig,Vs_sig] = get_segregate_pdf(xa,xv,prior_mu,A_sig,V_sig,prior_sig,eval_range)

if length(xa) == 1
        As_mu = (xa/A_sig^2 + prior_mu/prior_sig^2) / (1/A_sig^2 + 1/prior_sig^2); %incorporating prior sig on unisensory estimates
        Vs_mu = (xv/V_sig^2 + prior_mu/prior_sig^2) / (1/V_sig^2 + 1/prior_sig^2);
        As_sig = sqrt((1/A_sig^2 + 1/prior_sig^2)^-1);
        Vs_sig = sqrt((1/V_sig^2 + 1/prior_sig^2)^-1);
        A_pdf = normpdf(eval_range,As_mu,As_sig); %evaluate likelihood at each actual saccade location
        V_pdf = normpdf(eval_range,Vs_mu,Vs_sig);
        seg_pdf = (A_pdf+V_pdf)/2; %normalizing so that total probability sums to 1
else
        As_mu = (xa./A_sig^2 + prior_mu/prior_sig^2) / (1/A_sig^2 + 1/prior_sig^2); %incorporating prior sig on unisensory estimates
        Vs_mu = (xv./V_sig^2 + prior_mu/prior_sig^2) / (1/V_sig^2 + 1/prior_sig^2);
        As_sig = sqrt((1/A_sig^2 + 1/prior_sig^2)^-1);
        Vs_sig = sqrt((1/V_sig^2 + 1/prior_sig^2)^-1);
        expanded_range(1,1,1,:) = eval_range; %putting eval range in 4th dim to allow for integration later on.
        A_pdf = bsxfun_normpdf(expanded_range,As_mu,As_sig); %evaluate likelihood at each actual saccade location
        V_pdf = bsxfun_normpdf(expanded_range,Vs_mu,Vs_sig);
        seg_pdf = bsxfun(@plus,A_pdf,V_pdf)/2; %normalizing so that total probability sums to 1

end


end