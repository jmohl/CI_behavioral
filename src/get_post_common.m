%% calculate probability of common cause
%
% -------------------
% Jeff Mohl
% 6/22/18
% -------------------
%
% Description: returns  scalar probability of common cause between 0 and 1
%
% Inputs:
% A_mu, V_mu, prior_mu: means of A, V, and prior components
% A_sig, V_sig, prior_sig: std  of same


function post_common = get_post_common(A_mu,V_mu,prior_mu,A_sig,V_sig,prior_sig,prior_common)
%calculating post common
Pav_c1 = 1/(2*pi*sqrt(V_sig^2*A_sig^2 + V_sig^2*prior_sig^2 + prior_sig^2*A_sig^2))* ...
    exp(-0.5 * ((V_mu-A_mu)^2*prior_sig^2 + (V_mu-prior_mu)^2 * A_sig^2 + (A_mu-prior_mu)^2 * V_sig^2) / ...
    (V_sig^2*A_sig^2 + V_sig^2*prior_sig^2 + prior_sig^2*A_sig^2)); %eq 4 - A, V given one cause

Pav_c2 = 1/(2*pi*sqrt((V_sig^2 + prior_sig^2)*(prior_sig^2+A_sig^2)))* ...
    exp(-0.5 * ((V_mu-prior_mu)^2/(V_sig^2+prior_sig^2) + (A_mu-prior_mu)^2/(prior_sig^2+A_sig^2))); %eq 6 p 2 cause given av

post_common = Pav_c1 * prior_common /(Pav_c1 * prior_common + Pav_c2 * (1-prior_common)); %eq 2 - posterior on common cause
end