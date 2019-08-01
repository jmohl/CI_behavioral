%% plot_unity_combined
%
% -------------------
% Jeff Mohl
% 7/18/19
% -------------------
%
% Description: makes localization plots for specific example conditions.
% Will make subplots of a single figure for each condition in
% example_conds. meant to run on results from one subject, if want to run
% on multiple subjects collapsed into same plot will need to do something
% different
% 
function plot_localization(m,model,example_conds,plot_pred)
if length(m) > 1
   % pool data across model structures if provided with multiple 
   % otherwise, get data from single model structure
   % responses and fit distributions are 3d, so add another dimension to go
   % across subjects then take mean across this dimension.
   for subj_ind = 1:length(m)
       this_m = m{subj_ind};
       model_ind = ismember(vertcat(this_m.models{:}),model,'rows');
       conditions = this_m.conditions{model_ind};
       xrange = this_m.fitoptions.eval_midpoints;
       if model(3) == 3 %joint fit models have cell array rather than vectors for these
           responses(:,:,:,subj_ind) = this_m.responses{model_ind}{2};
           fit_dist(:,:,:,subj_ind) = this_m.fit_dist{model_ind}{2};           
       else
           responses(:,:,:,subj_ind) = this_m.responses{model_ind};
           fit_dist(:,:,:,subj_ind) = this_m.fit_dist{model_ind};
       end
       %get unimodal responses
       uni_ind = ismember(vertcat(this_m.models{:}),[0 0 4 1],'rows');
       uni_cond = this_m.conditions{uni_ind};
       A_resp(:,:,:, subj_ind) = this_m.responses{uni_ind}{1};
       V_resp(:,:,:, subj_ind) = this_m.responses{uni_ind}{2};
   end
   %get mean model and response fits
    responses = mean(responses,4);
    fit_dist = mean(fit_dist,4);
    A_resp = mean(A_resp,4);
    V_resp = mean(V_resp,4);
else
    %otherwise, get data from single model structure
    model_ind = ismember(vertcat(m.models{:}),model,'rows');
    conditions = m.conditions{model_ind};
    xrange = m.fitoptions.eval_midpoints;
    if model(3) == 3 %joint fit models have cell array rather than vectors for these
        responses = m.responses{model_ind}{2};
        fit_dist = m.fit_dist{model_ind}{2};
    else
        responses = m.responses{model_ind};
        fit_dist = m.fit_dist{model_ind};
    end
    %get unimodal responses
    uni_ind = ismember(vertcat(m.models{:}),[0 0 4 1],'rows');
    uni_cond = m.conditions{uni_ind};
    A_resp = m.responses{uni_ind}{1};
    V_resp = m.responses{uni_ind}{2};
end

%make figure
figure
plot_ind = 1;
for this_ex = example_conds
    subplot(1,length(example_conds),plot_ind)
    hold on
     plot_modelhist(responses(this_ex,:,:),fit_dist(this_ex,:,:),xrange,plot_pred);
    % add unimodal conditions to plot
    norm_A = A_resp(uni_cond{1}(:) == conditions(this_ex,1),:);
    norm_A = norm_A/sum(norm_A);%normalizing to probability, *2 because double saccade trials are rescaled this way for comparison with single saccade
    norm_V = V_resp(uni_cond{2}(:) == conditions(this_ex,2),:);
    norm_V = norm_V/sum(norm_V);%normalizing to probability,
    % get mean values
    mean_A = sum(norm_A.*xrange,2)./sum(norm_A,2);
    
    plot(xrange,norm_A,'--','Color',[1 0 0 .5]);
    plot(xrange,norm_V,'--','Color',[0 0 1 .5]);
%     text(mean_A,max(norm_A),sprintf('|%2.2f',mean_A))
if length(m) == 1
    title(sprintf('A=%d, V=%d, %s',conditions(this_ex,1:2),m.subject));
else
    title(sprintf('A=%d, V=%d',conditions(this_ex,1:2)));
end
     xlim([-25 25])
    plot_ind = 1 + plot_ind;
end
if plot_pred
    legend('Single sac, C=1','A sac, C=2','V sac, C=2','Modeled', 'A Unimodal', 'V Unimodal','Location','Best');
else
    legend('Single sac, C=1','A sac, C=2','V sac, C=2', 'A Unimodal', 'V Unimodal','Location','Best');
end

end