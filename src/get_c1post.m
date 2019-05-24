%% calculate posterior probability of common cause analytically
%
% -------------------
% Jeff Mohl
% 2/14/19
% -------------------
%
% Description: if inputs xa and xv are scalar, returns scalar
% probability of common cause between 0 and 1. If mu inputs are vectors,
% returns an array of probabilities for every combination of xa, xv. This
% is for the analytically solvable case where the prior is a single normal
% distribution.
%
% Inputs:
% A_mu, V_mu, prior_mu: means of A, V, and prior components
% A_sig, V_sig, prior_sig: std  of same
% Outputs:
% c1post: posterior probability (distribution) of common cause for xa and xv


function c1post = get_c1post(xa,xv,prior_mu,A_sig,V_sig,prior_sig,p_common)
%calculating post common
if length(xa) == 1
    Pav_c1 = 1/(2*pi*sqrt(V_sig^2*A_sig^2 + V_sig^2*prior_sig^2 + prior_sig^2*A_sig^2))* ...
        exp(-0.5 * ((xv-xa)^2*prior_sig^2 + (xv-prior_mu)^2 * A_sig^2 + (xa-prior_mu)^2 * V_sig^2) / ...
        (V_sig^2*A_sig^2 + V_sig^2*prior_sig^2 + prior_sig^2*A_sig^2)); %eq 4 - A, V given one cause
    
    Pav_c2 = 1/(2*pi*sqrt((V_sig^2 + prior_sig^2)*(prior_sig^2+A_sig^2)))* ...
        exp(-0.5 * ((xv-prior_mu)^2/(V_sig^2+prior_sig^2) + (xa-prior_mu)^2/(prior_sig^2+A_sig^2))); %eq 6 p 2 cause given av
    
    c1post = Pav_c1 * p_common /(Pav_c1 * p_common + Pav_c2 * (1-p_common)); %eq 2 - posterior on common cause
else
    Pav_c1 = 1/(2*pi*sqrt(V_sig^2*A_sig^2 + V_sig^2*prior_sig^2 + prior_sig^2*A_sig^2))* ...
        exp(-0.5 * (bsxfun(@minus, xv,xa).^2*prior_sig^2 + bsxfun(@plus,(xv-prior_mu).^2 * A_sig^2, (xa-prior_mu).^2 * V_sig^2)) / ...
        (V_sig^2*A_sig^2 + V_sig^2*prior_sig^2 + prior_sig^2*A_sig^2));%eq 4 - A, V given one cause
    
    Pav_c2 = 1/(2*pi*sqrt((V_sig^2 + prior_sig^2)*(prior_sig^2+A_sig^2)))* ...
        exp(-0.5 * (bsxfun(@plus,(xv-prior_mu).^2/(V_sig^2+prior_sig^2), (xa-prior_mu).^2/(prior_sig^2+A_sig^2))));
    
    c1post = Pav_c1 .* p_common ./(Pav_c1 .* p_common + Pav_c2 .* (1-p_common)); %eq 2 - posterior on common cause
    
end
end