%% Collecting adapted code from acerbi
%
% -------------------
% Jeff Mohl
% 2/14/19
% -------------------
%
% Description: a lot of code I tried to adapt from acerbi and ended up
% discarding. Here I am collecting that code in case I need to come back to
% it.

% calculate the C=2 PDF at each value of xrange, joint probability of xa xv
%NOTE: this code is adapted from lines 357:371 in
%VestBMS_BimodalLeftRightDatalike.m by luigi acerbi
%NOTE2: bsxfun_normpdf and _normcdf are also written by acerbi, but I have
%copied them into this project src
muc2_V = xrange_V.*prior_sig^2/(V_sig^2 + prior_sig^2);
sigmac2_V = V_sig*prior_sig/sqrt(V_sig^2 + prior_sig^2);
muc2_A = xrange_A.*prior_sig^2/(A_sig^2 + prior_sig^2);
sigmac2_A = A_sig*prior_sig/sqrt(A_sig^2 + prior_sig^2);
int_V = (bsxfun_normcdf(MAXRNG,muc2_V,sigmac2_V) - bsxfun_normcdf(-MAXRNG,muc2_V,sigmac2_V)); %JM what is the purpose of this step?
int_A = (bsxfun_normcdf(MAXRNG,muc2_A,sigmac2_A) - bsxfun_normcdf(-MAXRNG,muc2_A,sigmac2_A));
int_V = int_V .* bsxfun_normpdf(xrange_V,0,sqrt(V_sig^2 + prior_sig^2));
int_A = int_A .* bsxfun_normpdf(xrange_A,0,sqrt(A_sig^2 + prior_sig^2));
likec2 = bsxfun(@times, int_V, int_A);

% calculate the C=1 PDF
%NOTE: this code is adapted from lines 399:410 in
%VestBMS_BimodalLeftRightDatalike.m by luigi acerbi
mutilde = bsxfun(@plus, xrange_A.*A_sig^2, xrange_V.*V_sig^2)./(A_sig^2 + V_sig^2);
sigma2tilde = A_sig^2.*V_sig^2./(A_sig^2 + V_sig^2);
mucdf = mutilde.*prior_sig^2./(sigma2tilde + prior_sig^2);
sigmacdf = sqrt(sigma2tilde./(sigma2tilde + prior_sig^2))*prior_sig;
intc1 = (bsxfun_normcdf(MAXRNG, mucdf, sigmacdf) - bsxfun_normcdf(-MAXRNG, mucdf, sigmacdf));
likec1 = intc1 .* bsxfun_normpdf(xrange_A,xrange_V,sqrt(A_sig^2 + V_sig^2)) .* ...
    bsxfun_normpdf(mutilde,0,sqrt(sigma2tilde +prior_sig^2));
% get posterior for (C=1|xv,xa)
c1post = (likec1 * p_common)./(likec1 * p_common + (1-p_common)*likec2);
