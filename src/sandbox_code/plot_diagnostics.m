%% diagnostic plots for individual subjects
%
% -------------------
% Jeff Mohl
% 8/2/18
% -------------------
%
% Description: collecting diagnostic plots, currently only plotting biases
% in auditory and visual localization for evaluation of individual subjects

%% Is there a pattern of biases for auditory and visual guided saccades?
% relying on data from the get_unimodal_est function, stored in the
% fixed_parameters structure.
%auditory
function plot_diagnostics(fixed_params,data,subject)

%get unimodal aud estimates, little hacky
A_mu = [];
for i = 1:length(fixed_params.A_tars)
    A_tar = fixed_params.A_tars(i);
    this_data = table2array(data(data.A_tar == A_tar & isnan(data.V_tar),'valid_endpoints'));
    this_data = cell2mat(this_data);
    A_mu(end+1) = mean(this_data(:,1));
end


figure()
subplot(2,1,1)
plot(fixed_params.A_tars,A_mu, 'k.')
title('Auditory localization bias')
hold on
plot([-24,24], [-24,24], 'Color', [.5 .5 .5 .5])
xlabel('True target location')
ylabel('Mean response location')
xlim([-25,25])
ylim([-25,25])
subplot(2,1,2)
error_vector = fixed_params.A_tars-A_mu';
bar(fixed_params.A_tars,error_vector)
title('Error by target location')
xlabel('Tar Loc')
ylabel('Degrees error')

% %visual
% figure()
% subplot(2,1,1)
% plot(fixed_params.V_tars,fixed_params.V_mu, 'k.')
% title('Visual localization bias')
% hold on
% plot([-35,35], [-35,35], 'Color', [.5 .5 .5 .5])
% xlabel('True target location')
% ylabel('Mean response location')
% xlim([-35,35])
% ylim([-35,35])
% subplot(2,1,2)
% error_vector = fixed_params.V_tars-fixed_params.V_mu';
% bar(fixed_params.V_tars,error_vector)
% title('Error by target location')
% xlabel('Tar Loc')
% ylabel('Degrees error')

%notes: visual bias is quite small and seems additive (everything is
%slightly shifted left, matches what is fit in the model that includes a bias term)

%% find best fit for this unisensory bias, expect sigmoidal
% 'sigmoidal aud fit', generated using the curve fitting toolbox.
[xData, yData] = prepareCurveData( fixed_params.A_tars, A_mu');

% Set up fittype and options.
ft = fittype( 'a/(1+exp(-b*x))+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.495055162094339 0.753484628547185 0.30168010595327];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'sigmoidal aud fit' );
h = plot( fitresult, xData, yData );
hold on;
plot([-24,24], [-24,24], 'Color', [.5 .5 .5 .5])
legend('A mu resp', 'sigmoidal aud fit', 'unity ref line', 'Location', 'NorthWest');
% Label axes
xlabel('A tars')
ylabel('A mus')
grid on
title(sprintf('Unimodal Auditory localization bias - %s',subject),'Interpreter','none');
saveas(gcf,sprintf('results\\aud_bias_fits\\%s',subject),'png');
end
%notes: from this its very clear that there is a sigmoidal + additive bias
%here. The R^2 is like .99 on juno data

