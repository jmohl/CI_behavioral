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

%notes taken from looking at these plots
% - yoko has better overall behavioral performance on singles, especially with
% localizing auditory targets (sig ~4, compressive bias ~8 deg on 24 deg
% targets)
% - yoko has many more lapse trials, which result in things like major
% outlier V saccades
% - Juno sig ~8, compressive bias and leftward bias (mean for +24 is 10, but there are a lot of clear outliers)
% - for juno the AV-A sac conditions actually look pretty good.
% - Outliers seem to have a disproportionate effect on the standard
% deviation, and dealing with these might improve my fits
% - for come conditions, the "single" saccade distributions are clearly
% bimodal, indicating a lapse rather than a true single report (big problem
% for yoko, +12 vis conditions)
% All this looks like 2 potential problems to fix
% 1. A sacs or V sacs that are very far from their respective target should
% not be included. Can either make this a hard cutoff or based on the
% actual response distributions
% 2. There are a lot of "single" saccade responses that are not at all
% reflective of fusion, but more like a failure to initiate the seccond
% saccade, especially for yoko.

%% set initial state
local_directory = 'C:\Users\jtm47\Documents\Projects\CI_behavioral\';
cd(local_directory)
addpath('src', 'src\plotting','data');

data = load('data\Yoko_combined.mat');
data=data.tidy_data;
data = accuracy_filter(data);
data = data(logical(data.valid_tr),:);
data = data(data.n_sacs > 0,:);

%get resposne and condition vectors
A_data = data(strcmp(data.trial_type, 'A'),:);
V_data = data(strcmp(data.trial_type, 'V'),:);
AV_data = data(strcmp(data.trial_type, 'AV'),:);

[A_g, A_tars] = findgroups(A_data.A_tar);
[V_g, V_tars] = findgroups(V_data.V_tar);


A_responses =[A_data.A_tar,vertcat(A_data.valid_endpoints{:})];
A_responses = A_responses(:,1:2); % only include x loc


V_responses =[V_data.V_tar,vertcat(V_data.valid_endpoints{:})];
V_responses = V_responses(:,1:2); % only include x loc

% break out the responses to A, V and combined saccades to fit with normals
AV_A_resp = [AV_data(AV_data.n_sacs > 1,:).A_tar,AV_data(AV_data.n_sacs > 1,:).V_tar,vertcat(AV_data.A_endpoints{:})];
AV_V_resp = [AV_data(AV_data.n_sacs > 1,:).A_tar,AV_data(AV_data.n_sacs > 1,:).V_tar,vertcat(AV_data.V_endpoints{:})];
AV_C_resp = [AV_data(AV_data.n_sacs == 1,:).A_tar,AV_data(AV_data.n_sacs == 1,:).V_tar,vertcat(AV_data(AV_data.n_sacs ==1,:).valid_endpoints{:})];

[AV_A_g, AV_A_tars] = findgroups(AV_A_resp(:,1));
[AV_V_g, AV_V_tars] = findgroups(AV_V_resp(:,2));
[AV_C_g, C_tars] = findgroups(AV_C_resp(:,2));


tempfit = @(x)fitdist(x,'normal');
A_fits = splitapply(tempfit, A_responses(:,2),A_g);
V_fits = splitapply(tempfit, V_responses(:,2),V_g);
AV_A_fits = splitapply(tempfit, AV_A_resp(:,3),AV_A_g);
AV_V_fits = splitapply(tempfit, AV_V_resp(:,3),AV_V_g);
AV_C_fits = splitapply(tempfit, AV_C_resp(:,3),AV_C_g);

%find outliers in V data
% temp_outlier = @(x1){isoutlier(x1)};
% test3 = splitapply(temp_outlier,AV_V_resp(:,3),AV_V_g);
% cellfun(@sum,test3)
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
plot([A_tars(this_fit) A_tars(this_fit)],[0 max(A_plot)],'--k','LineWidth',2)
title(sprintf('%d A \n %2.2f mu, %2.2f sig',A_tars(this_fit),A_fits(this_fit).ParameterValues))
legend('Norm fit','Saccades','target loc')
end

%% AV condition, split up by target
for this_fit = 1:length(AV_A_fits)
figure
A_plot = pdf(AV_A_fits(this_fit),plot_range);
plot(plot_range,A_plot,'k')
hold on
histogram(AV_A_fits(this_fit).InputData.data,plot_range,'Normalization','probability')
plot([AV_A_tars(this_fit) AV_A_tars(this_fit)],[0 max(A_plot)],'--k','LineWidth',2)
title(sprintf('%d AV A sac \n %2.2f mu, %2.2f sig',AV_A_tars(this_fit),AV_A_fits(this_fit).ParameterValues))
legend('Norm fit','Saccades','target loc')
end

for this_fit = 1:length(AV_V_fits)
figure
V_plot = pdf(AV_V_fits(this_fit),plot_range);
plot(plot_range,A_plot,'k')
hold on
histogram(AV_V_fits(this_fit).InputData.data,plot_range,'Normalization','probability')
plot([AV_V_tars(this_fit) AV_V_tars(this_fit)],[0 max(V_plot)],'--k','LineWidth',2)
title(sprintf('%d AV V sac \n %2.2f mu, %2.2f sig',AV_V_tars(this_fit),AV_V_fits(this_fit).ParameterValues))
legend('Norm fit','Saccades','target loc')
end

for this_fit = 1:length(AV_C_fits)
figure
A_plot = pdf(AV_C_fits(this_fit),plot_range);
plot(plot_range,A_plot,'k')
hold on
histogram(AV_C_fits(this_fit).InputData.data,plot_range,'Normalization','probability')
plot([AV_V_tars(this_fit) AV_V_tars(this_fit)],[0 max(A_plot)],'--k','LineWidth',2)
title(sprintf('%dV AV single \n %2.2f mu, %2.2f sig',AV_V_tars(this_fit),AV_C_fits(this_fit).ParameterValues))
legend('Norm fit','Saccades','target loc')
end