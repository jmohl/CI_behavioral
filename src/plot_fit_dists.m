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
% data - raw data for this subject
% fit params, including ...
% A_mu, V_mu, prior_mu: means of A, V, and prior components
% A_sig, V_sig, prior_sig: std  of same
% prior_common: prior percent common cause
% plot_points: specify points for drawing the pdf (x coords and how smooth)
%
%
function plot_fit_dists(data,AV_pair,model,fit_params,fixed_params,plot_points)

%% extract params and data
prior_mu = fixed_params.prior_mu;
%target locations
xa = AV_pair(1);
xv = AV_pair(2);
%fit parameters
V_sig = fit_params.V_sig; %vis target sigma
Ac_sig = fit_params.Ac_sig; %close aud target sigma
Af_sig = fit_params.Af_sig; %far aud target sigma
if strcmp(model,'CI')
    prior_sig = fit_params.prior_sig;%sigma of centrality prior
    prior_common = fit_params.prior_common;%prior on common cause
end
if (abs(xa) > 12)
    A_sig = Af_sig;         %switch which auditory sigma is used, based on target location
else
    A_sig = Ac_sig;
end

%subset to only relevant trials
AVdata = vertcat(data(data.A_tar == xa & data.V_tar == xv,:).valid_endpoints{:,1});
AVdata = AVdata(:,1);
Adata = vertcat(data(data.A_tar == xa & isnan(data.V_tar),:).valid_endpoints{:,1});
Adata = Adata(:,1);
Vdata = vertcat(data(data.V_tar == xv & isnan(data.A_tar),:).valid_endpoints{:,1});
Vdata = Vdata(:,1);


%% get pdfs

% get pdf assuming integration
int_pdf = get_integrate_pdf(xa,xv,prior_mu,A_sig,V_sig,prior_sig,plot_points);

%get pdf assuming segregation
[seg_pdf,A_pred,V_pred] = get_segregate_pdf(xa,xv,prior_mu,A_sig,V_sig,prior_sig,plot_points);

%get posterior on common cause
post_common = get_post_common(xa,xv,prior_mu,A_sig,V_sig,prior_sig,prior_common);

%combined based on CI pdf
CI_pdf = post_common*int_pdf + (1-post_common)*seg_pdf;

%% plotting
figure()
set(gcf,'Position',[100,60,1049,895])

%unimodal plots
subplot(4,1,1)
hold on
histogram(Adata,plot_points,'FaceColor','r','FaceAlpha',.25, 'Normalization','Probability')
plot(plot_points,A_pred,'r')
histogram(Vdata,plot_points,'FaceColor','b','FaceAlpha',.25, 'Normalization','Probability')
plot(plot_points,V_pred,'b')
title(sprintf('unimodal responses: %d A, %d V',xa,xv))
ylabel('p')

% integrate plots
subplot(4,1,2)
hold on
histogram(AVdata,plot_points,'FaceColor','k','FaceAlpha',.25, 'Normalization','Probability')
plot(plot_points,int_pdf,'k')
title('integrate only')
ylabel('p')

% segregate plots
subplot(4,1,3)
hold on
histogram(AVdata,plot_points,'FaceColor','k','FaceAlpha',.25, 'Normalization','Probability')
plot(plot_points,seg_pdf,'k')
title('segregate only')
ylabel('p')

% CI plot
subplot(4,1,4)
hold on
histogram(AVdata,plot_points,'FaceColor','k','FaceAlpha',.25, 'Normalization','Probability')
plot(plot_points,CI_pdf,'k')
text(24,.2,sprintf('post common: %0.2f',post_common));
title('Full model')
ylabel('p')
xlabel('position (deg)')

end