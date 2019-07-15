%% plot_unity_combined
%
% -------------------
% Jeff Mohl
% 4/24/19
% -------------------
%
% Description: similar to plot_unity, but combines across subjects and
% model fits
% 

function plot_unity_combined(model_array,model)

hold on
for ind = 1:length(model_array)
    m=model_array{ind};
    model_ind = ismember(vertcat(m.models{:}),model,'rows');
    conditions = m.conditions{model_ind};
    if model(3) == 3 %joint fit models have cell array rather than vectors for these
        responses = m.responses{model_ind}{1};
        fit_dist = m.fit_dist{model_ind}{1};
    else
        responses = m.responses{model_ind};
        fit_dist = m.fit_dist{model_ind};
    end
    deltaAV = conditions(:,1)-conditions(:,2);
    %convert responses to percentages
    responses = responses./ sum(responses,2);
    %group by delta values
    [gav,avlabels] = findgroups(deltaAV);
    mean_resp(:,ind) = splitapply(@mean, responses(:,1),gav);
    mean_fit(:,ind) = splitapply(@mean, fit_dist(:,1),gav);
    
    % plot individual subject lines
%     plot(deltaAV,responses(:,1),'.k','MarkerSize',5) %was plotting each condition as a dot, but don't really like it
    plot(avlabels,mean_resp(:,ind),'Color', [.75 .75 .75 .5],'LineWidth',.5);

end
%across subject means and SEM
subj_mean = mean(mean_resp,2);
subj_sem = std(mean_resp,1,2)/sqrt(length(model_array));
model_mean = mean(mean_fit,2);
l1 = errorbar(avlabels,subj_mean,subj_sem,'k','LineWidth',2);
l2 = plot(avlabels,model_mean,'g','LineWidth',2);

hold off
legend([l1 l2],'Mean Response','Model predicted')
xlabel('\Delta AV')
ylabel('% unity judgement')
title('% of trials reported unity by target separation')

end