%% get unisensory prior
%
% -------------------
% Jeff Mohl
% 5/24/19
% -------------------
%
% Description: get the unisensory prior under the given model assumptions,
% for both A and V conditions
%
% Inputs:
% prior_type: 
%   2: discrete uniform prior over locations
%   3: balanced gaussian mixture prior with means(locations) and with std sigma_prior
% Outputs:
% prior: 1x(xrange) probability mass matrix 

function prior = get_unisensory_prior(prior_type,xrange,locations,sigma)
switch prior_type
    case 1 %naive normal prior
        prior = normpdf(xrange,locations,sigma);
    case 2 %uniform prior, over visual target locations
        prior = zeros(size(xrange));
        prior(ismember(xrange,locations)) = 1/length(locations);
    case 3 %mixture of gaussians prior
        prior = zeros(size(xrange));
        for loc = locations
            prior = prior+normpdf(xrange,loc,sigma)/length(locations);
        end             
end 
prior = prior/sum(prior,'all'); %normalize so that prior sums to 1, to deal with potentially over-broad priors (should treat them as uniform over range at the limit)
end
