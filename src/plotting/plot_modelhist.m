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

function plot_modelhist(saccades,predicted,xlocs)
figure;

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
    norm_predicted = predicted*abs(xlocs(1)-xlocs(2))^2;
    %single saccades will be along the diagonal
    I_mat = logical(eye([length(xlocs),length(xlocs)]));
    sing_sacs = norm_saccades(:,I_mat);
    norm_saccades(:,I_mat) = 0; %remove saccades that are AV from counts
    A_sacs = sum(norm_saccades,2)/2; %divide by 2 to make probabilities sum to 1
    V_sacs = sum(norm_saccades,3)/2;
    sac_bar = bar(xlocs,[sing_sacs(:),A_sacs(:),V_sacs(:)],'stacked');%normalizing to probability
    sac_bar(1).FaceColor = [.2 .2 .2];
    sac_bar(2).FaceColor = [1 .5 .5];
    sac_bar(3).FaceColor = [.5 .5 1];
    hold on
    projected_pred = (squeeze(sum(norm_predicted,2)) + squeeze(sum(norm_predicted,3))')/2; %divide by 2 to normalize
    plot(xlocs,projected_pred,'LineWidth',1.5,'Color',[.5 .5 .5]); %scaling by bin width so that the probability is matched correctly
    xlim([-40 40])
    legend('Single Sac','A sac','V sac','Modeled','Location','Best');
    xlabel('endpoint location (degrees)');
    ylabel('% saccades in bin')
    set(gca,'box','off')
    %probably need to double check this rsq calc cause I'm getting weird
    %values.
%     Rsq = 1-sum((norm_saccades(:)-norm_predicted(:)).^2)/sum((norm_saccades(:)-mean(norm_saccades(:))).^2); 
%     text(-39,.1,sprintf('R^2:%0.2f',Rsq));
end
    
end