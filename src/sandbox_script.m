%% sandbox script for CI project
%
% -------------------
% Jeff Mohl
% 7/30/18
% -------------------
%
% Description: for holding things that are in progress or just used to
% satisfy temporary curiosity

%% Is there a pattern of biases for auditory and visual guided saccades?
% relying on data from the get_unimodal_est function, stored in the
% fixed_parameters structure.
%auditory
figure()
subplot(2,1,1)
title('Auditory localization bias')
plot(fixed_params.A_tars,fixed_params.A_mu, 'k.')
hold on
plot([-24,24], [-24,24], 'Color', [.5 .5 .5 .5])
xlabel('True target location')
ylabel('Mean response location')
xlim([-25,25])
ylim([-25,25])
subplot(2,1,2)
title('Error by target location')
error_vector = fixed_params.A_tars-fixed_params.A_mu';
bar(fixed_params.A_tars,error_vector)

%visual
figure()
subplot(2,1,1)
title('Visuallocalization bias')
plot(fixed_params.V_tars,fixed_params.V_mu, 'k.')
hold on
plot([-35,35], [-35,35], 'Color', [.5 .5 .5 .5])
xlabel('True target location')
ylabel('Mean response location')
xlim([-35,35])
ylim([-35,35])
subplot(2,1,2)
title('Error by target location')
error_vector = fixed_params.V_tars-fixed_params.V_mu';
bar(fixed_params.V_tars,error_vector)

%notes: visual bias is quite small and seems additive (everything is
%slightly shifted left, matches what is fit in the model that includes a bias term)

%% find best fit for this unisensory bias, expect sigmoidal
% 'sigmoidal aud fit', generated using the curve fitting toolbox.
[xData, yData] = prepareCurveData( fixed_params.A_tars, fixed_params.A_mu');

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
title('Unimodal Auditory localization bias - H02 (600tr)');
saveas(gcf,sprintf('results\\aud_bias_fit_H02'),'png');

%notes: from this its very clear that there is a sigmoidal + additive bias
%here. The R^2 is like .99. This was done on Juno data. I'm going to try
%to do the same thing with human data and see how it looks.

