%% calculate probability of common cause
%
% -------------------
% Jeff Mohl
% 2/14/19
% -------------------
%
% Description: if inputs xa and xv are scalar, returns scalar
% probability of common cause between 0 and 1. If mu inputs are vectors,
% returns an array of probabilities for every combination of xa, xv.
%
% Inputs:
% A_mu, V_mu, prior_mu: means of A, V, and prior components
% A_sig, V_sig, prior_sig: std  of same


function post_common = get_post_common(xa,xv,prior_mu,A_sig,V_sig,prior_sig,prior_common)
%calculating post common
if length(xa) == 1
    Pav_c1 = 1/(2*pi*sqrt(V_sig^2*A_sig^2 + V_sig^2*prior_sig^2 + prior_sig^2*A_sig^2))* ...
        exp(-0.5 * ((xv-xa)^2*prior_sig^2 + (xv-prior_mu)^2 * A_sig^2 + (xa-prior_mu)^2 * V_sig^2) / ...
        (V_sig^2*A_sig^2 + V_sig^2*prior_sig^2 + prior_sig^2*A_sig^2)); %eq 4 - A, V given one cause

    Pav_c2 = 1/(2*pi*sqrt((V_sig^2 + prior_sig^2)*(prior_sig^2+A_sig^2)))* ...
        exp(-0.5 * ((xv-prior_mu)^2/(V_sig^2+prior_sig^2) + (xa-prior_mu)^2/(prior_sig^2+A_sig^2))); %eq 6 p 2 cause given av

    post_common = Pav_c1 * prior_common /(Pav_c1 * prior_common + Pav_c2 * (1-prior_common)); %eq 2 - posterior on common cause
else
    
    
end
end