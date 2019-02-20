%% correcting for eye tracker induced bias
%
% -------------------
% Jeff Mohl
% 8/2/18
% -------------------
%
% Description: eye tracker is calibrated manually during recording, but
% this process is not perfect. Thsi function effectively adjusts the gain
% and offset so that saccades are as well aligned with real space as
% possible while maintaining a linear relationship (taking visually guided
% saccades as approximately unbiased). It acts only on the x (horizontal)
% component of the data


function [data] = get_bias_corrected_data(data)

V_data = data(strcmp(data.trial_type,'V'),:);
%need to strip out occasional second saccades
saccades_vector = vertcat(V_data.valid_endpoints{:,1});
saccades_vector = saccades_vector(:,1);

%find linear coefficients for visual saccades
coeffs = polyfit(V_data.V_tar(:,1),saccades_vector,1);

for i = 1:height(data)
    this_endpoints = data.valid_endpoints{i};
    this_endpoints(:,1) = this_endpoints(:,1)/coeffs(1) - coeffs(2); %only adjusting first column, x component
    data.valid_endpoints{i}= this_endpoints;
    if data.n_sacs(i) > 1 && ismember('A_endpoints',data.Properties.VariableNames) %adjusted A and V saccades, if present
        this_A = data.A_endpoints{i};
        this_A(:,1) = this_A(:,1)/coeffs(1) - coeffs(2);
        this_V = data.V_endpoints{i};
        this_V(:,1) = this_V(:,1)/coeffs(1) - coeffs(2);
        data.A_endpoints{i} = this_A;
        data.V_endpoints{i} = this_V;
    end

end

end