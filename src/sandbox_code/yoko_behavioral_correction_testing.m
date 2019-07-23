%% yoko behavioral filtration 
%
% -------------------
% Jeff Mohl
% 7/23/19
% -------------------
%
% Description: Yoko has two behavioral deficiencies that are pretty
% annoying from a modeling perspective. 
% 1. She very often makes two saccades for some of the conditions which are
% actually coincident. These are mostly on the order of ~5 degrees. I might
% not get rid of it, but it sure is annoying. One solution would be to
% label saccades only based on the horizontal component (way back in
% tidy_data generation), but I don't really want to do that. 
% 2. She will often make a single saccade to either the aud or vis
% location, then give up. These should not really be interpreted as
% contributing to the "fused" distribution, but there is nothing in the
% model that prevents this. I have a "strict" behavior filter that corrects
% for this but, but it makes the unity judgement plot suspect since it's
% basically impossible to have a single saccade trial when the targets are
% well separated, so the line goes to 0 with 0 variance.


AV_data = tidy_data(strcmp(tidy_data.trial_type,'AV'),:);
% AV_data_strict = raw_data(strcmp(tidy_data.trial_type,'AV'),:);

% AV_data = AV_data_strict;
AV_data = AV_data(logical(AV_data.valid_tr),:);


AV_66 = AV_data(AV_data.A_tar == -6 & AV_data.V_tar == -6,:);

A_ep = vertcat(AV_66.A_endpoints{:});
V_ep = vertcat(AV_66.V_endpoints{:});
A_ep = A_ep(:,1);
V_ep = V_ep(:,1);

figure
plot(A_ep,V_ep,'k.')
xlabel('A endpoint')
ylabel('V endpoint')

figure
histogram(V_ep-A_ep)