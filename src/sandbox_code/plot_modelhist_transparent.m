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
% NOTE: 8/14/19 I spent some time using histograms instead of stacked bar
% plots, for the purposes of making things more clear. Unfortunately for
% directly comparing with the model this doesn't work, because the model is
% going to collapse across all the conditions so the probabilities won't
% match. If I want to do it using real pdfs then it has to be split between
% single saccade and dual saccade conditions. 

function plot_modelhist_transparent(saccades,predicted,xlocs, plot_pred)
global model_color aud_color vis_color

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
%     Rsq = 1-sum((norm_saccades-norm_predicted).^2)/sum((norm_saccades-mean(norm_saccades)).^2);
%     text(-39,.1,sprintf('R^2:%0.2f',Rsq));
elseif ndims(saccades) == 3

%convert back to saccade values, so can use histogram instead
    sac_data = saccades;
    binedges = xlocs(1:end-1)+.5;
    I_mat = logical(eye([length(xlocs),length(xlocs)]));
    sing_sacs = repelem(xlocs,squeeze(sac_data(:,I_mat)));
    sac_data(:,I_mat) = 0; %remove all single saccades from dist
    A_sacs = repelem(xlocs,squeeze(sum(sac_data,2)));
    V_sacs = repelem(xlocs,squeeze(sum(sac_data,3)));
    histogram(A_sacs,binedges,'FaceColor',aud_color)
    hold on
    histogram(V_sacs,binedges,'FaceColor',vis_color)
    histogram(sing_sacs,binedges,'FaceColor',[.1 .1 .1], 'FaceAlpha',.35)
    
    if plot_pred
        norm_predicted = predicted*abs(xlocs(1)-xlocs(2))^2;
        pred_sing = norm_predicted(:,I_mat) * (length(sing_sacs)+length(A_sacs)+ length(V_sacs));
        norm_predicted(:,I_mat) = 0;
        pred_A = squeeze(sum(norm_predicted,2))*(length(sing_sacs)+length(A_sacs)+ length(V_sacs));
        pred_V = squeeze(sum(norm_predicted,3))*(length(sing_sacs)+length(A_sacs)+ length(V_sacs));
        plot(xlocs,pred_sing,'--','LineWidth',1.5,'Color',[.25 .25 .25]); %scaling by bin width so that the probability is matched correctly
        plot(xlocs,pred_A,'--','LineWidth',1.5,'Color',aud_color); %scaling by bin width so that the probability is matched correctly
        plot(xlocs,pred_V,'--','LineWidth',1.5,'Color',vis_color); %scaling by bin width so that the probability is matched correctly

    end
    xlim([-40 40])
    xlabel('endpoint location (degrees)');
    ylabel('n saccades in bin')
    set(gca,'box','off')
end

end