%% plot the distance between A and V saccades as a function of target separation
%
% -------------------
% Jeff Mohl
% 4/24/19
% -------------------
%
% Description: plots for the localization of targets, but grouping by
% target separation and using the A to V separation as a metric. This is
% something that compares more directly with the plot_unity plots, but I'm
% not 100% sure I'll use it yet.

function plot_sac_by_sep(conditions,responses,fit_dist)

%using the data from the model fit structures, I am a little hampered
%because the data is binned in 1 degree bins rather than being actually
%continuous. If I want to not have this issue I need to go back to the raw
%data and get everything that way. But that doesn't solve the issue of what
%to do with the model fits, which are still going to be in the matrix
%format.

% This is a somewhat annoying problem because the separation between the
% two saccades is the difference between the A and V index, and you need to
% count up saccades along these non-diagonal diagonals.

delta_counts = zeros(size(conditions,1),101);
fit_counts = zeros(size(conditions,1),101);
for ic = 1:size(conditions,1) %condition number
    for delta = 0:100
        include_array = zeros(101,101);
        for k=1:101-delta
            j=k+delta;
            include_array(j,k) = 1;
        end
        include_array = logical(include_array + include_array');
        delta_counts(ic,delta+1) = sum(responses(ic,include_array))/sum(responses(ic,:,:),'all');
        fit_counts(ic,delta+1) = sum(fit_dist(ic,include_array),'all');
    end
end

deltaAV = abs(conditions(:,1)-conditions(:,2));

%convert delta counts to mean AV sep
all_sep = sum(delta_counts .* [0:100],2); %summing here because I have already normalized to probability
all_fit = sum(fit_counts .* [0:100],2);

[gav,gav_labs] = findgroups(deltaAV);
mean_sep = splitapply(@mean,all_sep,gav);
mean_fit = splitapply(@mean,all_fit,gav);

plot(deltaAV,all_sep,'k.')
hold on
plot(deltaAV,all_fit,'b.')
plot(gav_labs,mean_sep,'k','LineWidth',2);
plot(gav_labs,mean_fit,'b','LineWidth',2);
refline(1,0)
title('Mean distance between A and V saccades by distance between A and V targets')
xlabel('\Delta AV tars')
ylabel('\Delta AV saccades')
legend('Condition Responses','Condition fits','Mean Resp','Mean Fits','Location','NorthWest')

end