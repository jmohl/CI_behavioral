%% make box plot of saccades by target condition
% -------------------
% Jeff Mohl
% 1/03/19
% -------------------

% Description: make a blox plot showing the distribution of saccades for
% each target condition. Includes a reference line for actual target
% location. Meant only for single saccade trials

function plot_sac_box(tidy_data)
figure()
hold on
switch tidy_data.trial_type{1}
    case 'A'
        [g,glab] = findgroups(tidy_data.A_tar);
        title('Auditory Trials')
    case 'V'
        [g,glab] = findgroups(tidy_data.V_tar);
        title('Visual Trials')
    case 'AV'
%         if max(isnan(tidy_data.V_tar))
            error('sacbox is only meant for AV trials with A tar == V Tar')
%         else
%             [g,glab] = findgroups(tidy_data.V_tar);
%             title('AV same Trials')
%         end
end
%sac_vec = get_first_saccade(tidy_data);
sac_vec = vertcat(tidy_data.valid_endpoints{:});
boxplot(sac_vec(:,1),g,'Labels',glab); %only running on horizontal component of saccade)
ylim([-30,30]);
xlabel('Tar Location (deg)')
ylabel('Saccade Endpoints (deg)')
plot([glab,glab],'Color',[.5 .5 .5 .5]) %reference line
hold off

end