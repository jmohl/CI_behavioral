%% plot modeled vs actual saccade distributions
%
% -------------------
% Jeff Mohl
% 2/14/19
% -------------------
%
% Description: plot the modeled and response histograms for the saccades in
% a given condition. Haven't figured out a good way to collapse across
% conditions yet, maybe fit mixture of gaussians and compute r^2?
%
%

function plot_modelhist(saccades,predicted,xlocs, plot_pred)
global model_color aud_color vis_color

show_error = 0; %defaults to false, if predicted is a cell array, will assume true

if ndims(saccades) == 2 %format for location estimation when all saccades are put together
    norm_saccades=saccades/sum(saccades);
    norm_predicted = predicted*abs(xlocs(1)-xlocs(2));
    sac_bar = bar(xlocs,norm_saccades);%normalizing to probability
    sac_bar.FaceAlpha = 0.5;
    sac_bar.FaceColor = [125/255 39/255 125/255];
    hold on
    plot(xlocs,norm_predicted,'LineWidth',2,'Color','k'); %scaling by bin width so that the probability is matched correctly
    xlim([-40 40])
    legend('Actual','Modeled');
    xlabel('endpoint location (degrees)');
    ylabel('% saccades in bin')
    
    %adding R_squared to plots
    Rsq = 1-sum((norm_saccades-norm_predicted).^2)/sum((norm_saccades-mean(norm_saccades)).^2);
    text(-39,.1,sprintf('R^2:%0.2f',Rsq));
elseif ndims(saccades) == 3
    % this condition is a little trickier because it's 3 dimensional for
    % the analysis but I want to compress it down to 2 dimensions for
    % plotting. The predicted (likelihood) is simple because I can just
    % marginalize across A and V and sum those, but for the counts it
    % doesn't seem so simple. Maybe i can do the same thing after
    % normalizing the counts by probability?
    norm_saccades=saccades/sum(saccades(:));
    %single saccades will be along the diagonal
    I_mat = logical(eye([length(xlocs),length(xlocs)]));
    sing_sacs = norm_saccades(:,I_mat);
    norm_saccades(:,I_mat) = 0; %remove saccades that are AV from counts
    A_sacs = sum(norm_saccades,2); %divide by 2 to make probabilities sum to 1
    V_sacs = sum(norm_saccades,3); %JTM 7/18/19 made change here, removed division by 2
    sac_bar = bar(xlocs,[A_sacs(:),V_sacs(:),sing_sacs(:)],'stacked');%normalizing to probability
    sac_bar(1).FaceColor = aud_color;
    sac_bar(2).FaceColor = vis_color;
    sac_bar(3).FaceColor = [.2 .2 .2];
    hold on
    if plot_pred
        if length(predicted) == 2
            mean_predicted = predicted{1};
            sem_predicted = predicted{2};
            show_error = 1;
        else
            mean_predicted = predicted;
        end
        norm_predicted = mean_predicted*abs(xlocs(1)-xlocs(2))^2;
        norm_pred_sing = norm_predicted(:,I_mat);
        norm_predicted(:,I_mat) = 0;
        pred_A = squeeze(sum(norm_predicted,2));
        pred_V = squeeze(sum(norm_predicted,3));
        projected_pred = pred_A + pred_V' + norm_pred_sing';
        if show_error
            norm_sem = sem_predicted*abs(xlocs(1)-xlocs(2))^2;
            norm_sem_sing = norm_sem(:,I_mat);
            norm_sem(:,I_mat) = 0;
            norm_sem_A = squeeze(sum(norm_sem,2));
            norm_sem_V = squeeze(sum(norm_sem,3));
            norm_sem_pred = norm_sem_A + norm_sem_V' + norm_sem_sing';
            model_sem_bnd =[projected_pred + norm_sem_pred,projected_pred - norm_sem_pred];
            % plot confidence intervals for model
            x_plot =[xlocs, fliplr(xlocs)];
            y_plot=[model_sem_bnd(:,1); flipud(model_sem_bnd(:,2))];
            fill(x_plot, y_plot', 1,'facecolor', model_color, 'edgecolor', model_color, 'facealpha', 0.5);
        end
         plot(xlocs,projected_pred,'LineWidth',1.5,'Color',model_color); %scaling by bin width so that the probability is matched correctly
    end
    xlim([-40 40])
    xlabel('endpoint location (degrees)');
    ylabel('% saccades in bin')
    set(gca,'box','off')
    
end

end