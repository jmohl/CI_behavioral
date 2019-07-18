%% plot_unity
%
% -------------------
% Jeff Mohl
% 4/24/19
% -------------------
%
% Description: plots for the unity judgement fits, given a fit model
% structure and the desired model to plot.
% 
% For these plots, each condition is grouped with other conditions with the
% same A-V target separation, with the standard deviation and mean
% calculated within groups. When multiple subjects are combined, the
% conditions are grouped across subjects (rather than determining the
% mean/sd etc for each subject and then combining, which might be the better way to do it).

function plot_unity(conditions,responses,fit_dist)

%convert from target values in conditions to scalar separation value
deltaAV = conditions(:,1)-conditions(:,2);
%convert responses to percentages
responses = responses./ sum(responses,2);
%group by delta values
[gav,avlabels] = findgroups(deltaAV);
mean_resp = splitapply(@mean, responses(:,1),gav);
mean_fit = splitapply(@mean, fit_dist(:,1),gav);
mean_sd = splitapply(@std,responses(:,1),gav);

plot(deltaAV,responses(:,1),'.k','MarkerSize',5)
hold on
errorbar(avlabels,mean_resp,mean_sd,'k','LineWidth',2);
plot(avlabels,mean_fit,'g','LineWidth',2);

hold off
legend('Single Conditions','Mean Response','Model predicted')
xlabel('Target Separation (A - V)')
ylabel('% unity judgement')
title('% of trials reported unity by target separation')
end