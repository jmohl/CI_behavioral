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
global aud_color vis_color
plot_unimodal = 0;
plot_error = 1;    

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
    responses = sum(responses,4);
    fit_dist_sem = std(fit_dist,[],4)/sqrt(size(fit_dist,4));
    fit_dist = mean(fit_dist,4);
    A_resp = sum(A_resp,4);
    V_resp = sum(V_resp,4);
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
     plot_modelhist(responses(this_ex,:,:),{fit_dist(this_ex,:,:),fit_dist_sem(this_ex,:,:)},xrange,plot_pred);
         xlim([-10,30])

     %8/13/19 note: the way I'm doing all of this is somewhat janky. It might be
     %better to go back to the raw saccades rather than using the binned data,
     %but that sounds like a future issue.
     if plot_unimodal
         % add unimodal conditions to plot
         norm_A = A_resp(uni_cond{1}(:) == conditions(this_ex,1),:);
         norm_A = norm_A/sum(norm_A);%normalizing to probability, *2 because double saccade trials are rescaled this way for comparison with single saccade
         norm_V = V_resp(uni_cond{2}(:) == conditions(this_ex,2),:);
         norm_V = norm_V/sum(norm_V);%normalizing to probability,
         % get mean values
         mean_A = sum(norm_A.*xrange,2)./sum(norm_A,2);
         mean_V = sum(norm_V.*xrange,2)./sum(norm_V,2);
         
         %get means from AV trials
%          norm_saccades=responses(this_ex,:,:)/sum(responses(this_ex,:,:),'all');
%          %single saccades will be along the diagonal
%          I_mat = logical(eye([length(xrange),length(xrange)]));
%          sing_sacs = norm_saccades(:,I_mat);
%          norm_saccades(:,I_mat) = 0; %remove saccades that are AV from counts
%          A_sacs = squeeze(sum(norm_saccades,2));
%          V_sacs = squeeze(sum(norm_saccades,3));
%          
%          AV_mean_A = sum(A_sacs'.*xrange)/sum(A_sacs);
%          AV_mean_V = sum(V_sacs.*xrange)/sum(V_sacs);
%          AV_mean_sing = sum(sing_sacs.*xrange)/sum(sing_sacs);
%          % plot mean indicator - triangle
%          max_p = .35;
%          plot(mean_A,max_p,'v','MarkerSize',4, 'MarkerEdgeColor',aud_color)
%          plot(mean_V,max_p,'v','MarkerSize',4, 'MarkerEdgeColor',vis_color)
%          plot(AV_mean_A,max_p,'v','MarkerSize',4,'MarkerEdgeColor',aud_color,'MarkerFaceColor',aud_color)
%          plot(AV_mean_V,max_p,'v','MarkerSize',4,'MarkerEdgeColor',vis_color,'MarkerFaceColor',vis_color)
%          plot(AV_mean_sing,max_p,'v','MarkerSize',4,'MarkerEdgeColor',[.2 .2 .2],'MarkerFaceColor',[.2 .2 .2])
%          
         plot(xrange,norm_A,'--','Color',aud_color);
         plot(xrange,norm_V,'--','Color',vis_color);
     end
     %add target location plots
    plot([conditions(this_ex,1) conditions(this_ex,1)],[0 .35],'--','Color',aud_color);
    plot([conditions(this_ex,2) conditions(this_ex,2)],[0 .35],'--','Color',vis_color);
if length(m) == 1
    title(sprintf('A=%d, V=%d, %s',conditions(this_ex,1:2),m.subject));
else
    title(sprintf('A=%d, V=%d',conditions(this_ex,1:2)));
end
    plot_ind = 1 + plot_ind;
end
if plot_pred
    legend('A sac, C=2','V sac, C=2','Single sac, C=1','Model SEM','Model', 'A target', 'V target','Location','Best');
else
    legend('A sac, C=2','V sac, C=2','Single sac, C=1', 'A target', 'V target','Location','Best');
end

end