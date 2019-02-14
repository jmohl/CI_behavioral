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
sac_bar = bar(xlocs,saccades/sum(saccades));%normalizing to probability
sac_bar.FaceAlpha = 0.25;
hold on
plot(xlocs,predicted*abs(xlocs(1)-xlocs(2)),'LineWidth',2,'k'); %scaling by bin width so that the probability is matched correctly
xlim([-40 40])
legend('Actual','Modeled');
xlabel('endpoint location (degrees)');
ylabel('% saccades in bin')

end