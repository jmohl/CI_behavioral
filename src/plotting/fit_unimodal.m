%% fit unimodal trial
%
% -------------------
% Jeff Mohl
% 2/27/19
% -------------------
%
% Description: Code for fitting normal distributions to unimodal trials,
% for comparison with fits from the combined condition trials

%note: in the future I should update this to include the prior, as it is in
%the CI models, but for now I am just interested in getting some
%approximations of unisensory perception

%% set initial state
local_directory = 'C:\Users\jtm47\Documents\GitHub\CI_behavioral\';
cd(local_directory)
addpath('src', 'src\plotting','data');

data = load('data\Yoko_combined.mat');
data=data.tidy_data;
data = data(logical(data.valid_tr),:);
data = data(data.n_sacs > 0,:);

%get resposne and condition vectors
A_data = data(strcmp(data.trial_type, 'A'),:);
V_data = data(strcmp(data.trial_type, 'V'),:);

[A_g, A_tars] = findgroups(A_data.A_tar);
[V_g, V_tars] = findgroups(V_data.V_tar);

A_responses =[A_data.A_tar,vertcat(A_data.valid_endpoints{:})];
A_responses = A_responses(:,1:2); % only include x loc


V_responses =[V_data.V_tar,vertcat(V_data.valid_endpoints{:})];
V_responses = V_responses(:,1:2); % only include x loc

tempfit = @(x)fitdist(x,'normal');
A_fits = splitapply(tempfit, A_responses(:,2),A_g);
V_fits = splitapply(tempfit, V_responses(:,2),V_g);

%% make plots
%notes: outliers make my norm_pdf fits worse than expected, Yoko behavior
%looks way better than Juno on this, much more accurate.
plot_range = -40:40;
for this_fit = 1:length(V_fits)
figure
V_plot = pdf(V_fits(this_fit),plot_range);
plot(plot_range,V_plot,'k')
hold on
histogram(V_fits(this_fit).InputData.data,plot_range,'Normalization','probability')
plot([V_tars(this_fit) V_tars(this_fit)],[0 max(V_plot)],'--k','LineWidth',2)
title(sprintf('%d V \n %2.2f mu, %2.2f sig',V_tars(this_fit),V_fits(this_fit).ParameterValues))
legend('Norm fit','Saccades','target loc')
end

for this_fit = 1:length(A_fits)
figure
A_plot = pdf(A_fits(this_fit),plot_range);
plot(plot_range,A_plot,'k')
hold on
histogram(A_fits(this_fit).InputData.data,plot_range,'Normalization','probability')
plot([A_tars(this_fit) A_tars(this_fit)],[0 max(V_plot)],'--k','LineWidth',2)
title(sprintf('%d A \n %2.2f mu, %2.2f sig',A_tars(this_fit),A_fits(this_fit).ParameterValues))
legend('Norm fit','Saccades','target loc')
end

