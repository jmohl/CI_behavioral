%% plot_unity
%
% -------------------
% Jeff Mohl
% 4/24/19
% -------------------
%
% Description: plots for the unity judgement fits, given a fit model
% structure and the desired model to plot.

function plot_unity(conditions,responses,fit_dist)

%convert from target values in conditions to scalar separation value
deltaAV = abs(conditions(:,1)-conditions(:,2));
%convert responses to percentages
responses = responses./ sum(responses,2);
%group by delta values
[gav,avlabels] = findgroups(deltaAV);
mean_resp = splitapply(@mean, responses(:,1),gav);
mean_fit = splitapply(@mean, fit_dist(:,1),gav);

plot(deltaAV,responses(:,1),'k.')
hold on
plot(avlabels,mean_resp,'k','LineWidth',2);
plot(avlabels,mean_fit,'g','LineWidth',2);
hold off
legend('Single Conditions','Mean Response','Model predicted')
xlabel('\Delta AV')
ylabel('% unity judgement')
title('% of trials reported unity by target separation')
end