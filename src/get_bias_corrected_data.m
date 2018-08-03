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
% saccades as approximately unbiased).


function [corrected_V_data, corrected_A_data, corrected_AV_data] = ...
    get_bias_corrected_data(V_endpoint_data,A_endpoint_data, AV_endpoint_data, V_tars)

%make put all V saccades in vector format
sac_vector = vertcat(V_endpoint_data{:});

%create matching target vector
tar_vector = [];
for i = 1:length(V_endpoint_data)
    tar_vector = vertcat(tar_vector,repmat(V_tars(i),length(V_endpoint_data{i}),1));
end

%find linear coefficients for visual saccades
coeffs = polyfit(tar_vector,sac_vector,1);

%make corrected cell arrays
corrected_V_data = V_endpoint_data;
corrected_A_data = A_endpoint_data;
corrected_AV_data = AV_endpoint_data;
for i = 1:length(V_endpoint_data)
    corrected_V_data{i} = V_endpoint_data{i}/coeffs(1) - coeffs(2); %dividing and subtracting by coefficient values will center data on unity line, removing linear bias
end
for i = 1:length(A_endpoint_data)
    corrected_A_data{i} = A_endpoint_data{i}/coeffs(1) - coeffs(2);
end
for i = 1:length(AV_endpoint_data)
    corrected_AV_data{i} = AV_endpoint_data{i}/coeffs(1) - coeffs(2);
end

end