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

end