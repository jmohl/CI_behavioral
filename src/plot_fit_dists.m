%% plot distributions based on fit parameters
%
% -------------------
% Jeff Mohl
% 6/25/18
% -------------------
%
% Description: script for visualizing the distributions based on the fit
% parameters for causal inference model. Currently only using for that
% model and not set up for simpler versions.
%
% Inputs:
% A_mu, V_mu, prior_mu: means of A, V, and prior components
% A_sig, V_sig, prior_sig: std  of same
% prior_common: prior percent common cause

function plot_fit_dists(data,A_tar,V_tar,fit_params,fixed_params,plot_points)

%% extract params and data
prior_mu = fixed_params.prior_mu;

V_sig = fit_params.V_sig; %vis target sigma
Ac_sig = fit_params.Ac_sig; %close aud target sigma
Af_sig = fit_params.Af_sig; %far aud target sigma
prior_sig = fit_params.prior_sig;%sigma of centrality prior
prior_common = fit_params.prior_common;%prior on common cause

A_data = data(strcmp(data.trial_type,'A'),:);
A_endpoints = get_response_endpoints(A_data(A_data.A_tar == A_tar,:),1,100); % require fix, 100ms buffer
A_endpoints = A_endpoints(:,1);

V_data = data(strcmp(data.trial_type,'V'),:);
V_endpoints= get_response_endpoints(V_data(V_data.V_tar == V_tar,:),1,100); 
V_endpoints = V_endpoints(:,1);

AV_data = data(strcmp(data.trial_type,'AV'),:);
AV_endpoints = get_response_endpoints(AV_data(AV_data.A_tar == A_tar & AV_data.V_tar == V_tar,:),1,100); % require fix, 100ms buffer
AV_endpoints = AV_endpoints(:,1);

if (abs(A_tar) > 12)
    A_sig = Af_sig;         %switch which auditory sigma is used, based on target location
else
    A_sig = Ac_sig;
end

%% get pdfs

% get pdf assuming integration
int_pdf = get_integrate_pdf(A_tar,V_tar,prior_mu,A_sig,V_sig,prior_sig,plot_points);

%get pdf assuming segregation
[seg_pdf,A_pred,V_pred] = get_segregate_pdf(A_tar,V_tar,prior_mu,A_sig,V_sig,prior_sig,plot_points);

%get posterior on common cause
post_common = get_post_common(A_tar,V_tar,prior_mu,A_sig,V_sig,prior_sig,prior_common);

%combined based on CI pdf
CI_pdf = post_common*int_pdf + (1-post_common)*seg_pdf;

%% plotting
figure()
set(gcf,'Position',[100,60,1049,895])

%unimodal plots
subplot(4,1,1)
hold on
histogram(A_endpoints,plot_points,'FaceColor','r','FaceAlpha',.25, 'Normalization','Probability')
plot(plot_points,A_pred,'r')
histogram(V_endpoints,plot_points,'FaceColor','b','FaceAlpha',.25, 'Normalization','Probability')
plot(plot_points,V_pred,'b')
title(sprintf('unimodal responses: %d A, %d V',A_tar,V_tar))
ylabel('p')

% integrate plots
subplot(4,1,2)
hold on
histogram(AV_endpoints,plot_points,'FaceColor','k','FaceAlpha',.25, 'Normalization','Probability')
plot(plot_points,int_pdf,'k')
title('integrate only')
ylabel('p')

% segregate plots
subplot(4,1,3)
hold on
histogram(AV_endpoints,plot_points,'FaceColor','k','FaceAlpha',.25, 'Normalization','Probability')
plot(plot_points,seg_pdf,'k')
title('segregate only')
ylabel('p')

% CI plot
subplot(4,1,4)
hold on
histogram(AV_endpoints,plot_points,'FaceColor','k','FaceAlpha',.25, 'Normalization','Probability')
plot(plot_points,CI_pdf,'k')
title('Full model')
ylabel('p')
xlabel('position (deg)')

end