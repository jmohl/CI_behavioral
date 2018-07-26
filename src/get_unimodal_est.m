%% estimate paramters for unimodal percepts
%
% -------------------
% Jeff Mohl
% 7/26/18
% -------------------
%
% Description: this function takes in the saccade endpoints and for the
% unimodal conditions (A and V) and 2 vetors, one for
% estimated mean and one for estimated std (sigma).

function [mu_vector, sig_vector] = get_unimodal_est(endpoint_data)

for i = 1:size(endpoint_data,2)
    this_data = endpoint_data{i};
    mu_vector(i) = mean(this_data);
    sig_vector(i) = std(this_data);
end