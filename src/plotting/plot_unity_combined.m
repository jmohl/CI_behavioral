%% plot_unity_combined
%
% -------------------
% Jeff Mohl
% 7/16/19
% -------------------
%
% Description: similar to plot_unity, but combines across subjects and
% model fits
% 

function plot_unity_combined(model_array,model)
global model_color 

 subj_ids = cell(length(model_array),1);
 subj_responses = zeros(length(model_array),20);

hold on
for ind = 1:length(model_array)
    m=model_array{ind};
    subj_ids{ind} = m.subject;
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
     subj_responses(ind,:) = responses(:,1);
    % plot individual subject lines
%     plot(deltaAV,responses(:,1),'.k','MarkerSize',5) %was plotting each condition as a dot, but don't really like it
    plot(avlabels,mean_resp(:,ind),'Color', [.75 .75 .75 .75],'LineWidth',.5);

end
%across subject means and SEM
subj_mean = mean(mean_resp,2);
subj_sem = std(mean_resp,1,2)/sqrt(length(model_array));
model_mean = mean(mean_fit,2);
model_sem = std(mean_fit,1,2)/sqrt(length(model_array));
model_sem_bnd =[model_mean + model_sem,model_mean - model_sem]; 

l1 = errorbar(avlabels,subj_mean,subj_sem,'k','LineWidth',2);
l2 = plot(avlabels,model_mean,'--','color', model_color,'LineWidth',2);

% plot confidence intervals for model
x_plot =[avlabels; flipud(avlabels)];
y_plot=[model_sem_bnd(:,1); flipud(model_sem_bnd(:,2))];
fill(x_plot, y_plot, 1,'facecolor', model_color, 'edgecolor', 'none', 'facealpha', 0.25);


hold off
ylim([0,1]);
legend([l1 l2],'Mean Response','Model predicted')
xlabel('Target Separation (A - V)')
ylabel('% unity judgement')
title('% of trials reported unity by target separation')

%% repeated measures anova
  subj_responses = subj_responses';
anova_table = array2table(subj_responses(:),'VariableNames',{'p_sing'});
table_ind = 1;
for ind = 1:length(model_array)
anova_table.subj(table_ind:table_ind + 19) = subj_ids(ind);
anova_table.tar_sep(table_ind:table_ind + 19) = deltaAV;
table_ind = table_ind+20;
end
[p, table, stats] = anovan(anova_table.p_sing,{anova_table.subj,anova_table.tar_sep},'random',1);

end