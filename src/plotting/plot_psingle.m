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
function plot_psingle(data,modelfit)

if nargin < 2
    modelfit = [];
end

% get condition vectors
conds = unique(data(:,2:3),'rows');
respbins = unique(data(:,4));

% get reponse for each condition (1 column single counts, 1 double
% counts )
responses = zeros(length(conds),length(respbins));
for ic = 1:length(responses)
    responses(ic,1) = sum(data(:,2) == conds(ic,1)&data(:,3) == conds(ic,2) & data(:,4) == 1); %count single saccade trials in cond ci
    responses(ic,2) = sum(data(:,2) == conds(ic,1)&data(:,3) == conds(ic,2) & data(:,4) == 2); %count double saccade trials in cond ci
end
p_sing = responses(:,1)./(sum(responses(:,:),2));

figure()
hold on
A_tars = unique(conds(:,1)); %grouping target disparities
for ii = 1:length(A_tars)
    this_tar = A_tars(ii);
    %disp = conds(conds(:,1) == this_tar,2) - this_tar; %target disparity
    disp = conds(conds(:,1) == this_tar,2); %using visual target location ,will almost certainly change in future but this provides maximum information
    p(ii) = plot(disp,p_sing(conds(:,1) == this_tar),'LineWidth',2);
    plot(disp,modelfit(conds(:,1) == this_tar),'LineWidth',2,'Color','k','LineStyle','--');
end

legend(p,cellstr(num2str(A_tars)),'Location','bestoutside')
title('Single saccade ratio by target separation')
ylabel('percent single saccade')
xlabel('Vis target location')

end