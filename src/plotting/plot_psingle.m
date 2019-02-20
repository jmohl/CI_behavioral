%% Plot percent single saccade by condition vs predictions by model
%
% -------------------
% Jeff Mohl
% 2/1/19
% -------------------
%
% Description: make plots for evaluating the percentage of single saccade
% trials in each condition. Also will plot the model fit results (if
% included)
%
% Inputs: data array used for model fitting, modelfit probabilities for
% each condition (1 row per condition).
%
% todo: this was just thrown together so I could evaluate model fits, and
% is not the final thing I want to use to make figures. it needs a ton of
% improvement. Also it somewhat duplicates functionality of plot_p_single
% in the AVD exploration project.
function plot_psingle(responses,conds,modelfit)

if nargin < 2
    modelfit = [];
end

figure()
hold on
A_tars = unique(conds(:,1)); %split into 24 and 
for ii = 1:length(A_tars)
    this_tar = A_tars(ii);
    %disp = conds(conds(:,1) == this_tar,2) - this_tar; %target disparity
    v_tar = conds(conds(:,1) == this_tar,2); %using visual target location ,will almost certainly change in future but this provides maximum information
    p(ii) = plot(v_tar,responses(conds(:,1) == this_tar,1)./sum(responses(conds(:,1) == this_tar,:),2),'LineWidth',2);
    plot(v_tar,modelfit(conds(:,1) == this_tar),'LineWidth',2,'Color',p(ii).Color,'LineStyle','--');
end

legend(p,cellstr(num2str(A_tars)),'Location','bestoutside')
title('Single saccade ratio by target separation')
ylabel('percent single saccade')
xlabel('Vis target location')

end