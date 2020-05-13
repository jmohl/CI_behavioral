%% determine 3+ saccade rate
%
% -------------------
% Jeff Mohl
% 4/21/20
% -------------------
%
% Description: determining ratio of trials with >2 saccades (only counting
% trials that would otherwise be considered valid)

yoko_data = load('data\Yoko_combined.mat');
yoko_data = yoko_data.tidy_data;

yoko_data = yoko_data(logical(yoko_data.valid_tr),:);
p_3_yoko = sum(yoko_data.n_sacs >= 3)/height(yoko_data);
p_2_yoko = sum(yoko_data.n_sacs == 2)/height(yoko_data);


juno_data = load('data\juno_combined.mat');
juno_data = juno_data.tidy_data;

juno_data = juno_data(logical(juno_data.valid_tr),:);
p_3_juno = sum(juno_data.n_sacs >= 3)/height(juno_data);
p_2_juno = sum(juno_data.n_sacs == 2)/height(juno_data);

h_data = [];
for ind = 2:8
    this_file = dir(sprintf('data\\H0%d_*',ind));
    this_h_data = load(this_file.name);
    this_h_data = this_h_data.tidy_data;
    h_data = [h_data;this_h_data];
end
        
h_data = h_data(logical(h_data.valid_tr),:);
p_3_h = sum(h_data.n_sacs >= 3)/height(h_data);
p_2_h = sum(h_data.n_sacs == 2)/height(h_data);


%plot third saccade locations, to determine if centered around fixation
%light or around target locations
h_data.spkdata = [];
juno_data.is_multiunit = [];
yoko_data.is_multiunit = [];

figure
hold on
trip_data = vertcat(h_data(h_data.n_sacs > 3,:),juno_data(juno_data.n_sacs > 3,:),yoko_data(yoko_data.n_sacs >3,:));
for ind = 1:height(trip_data)
    this_sac = trip_data(ind,:).valid_endpoints{1}(3,:);
scatter(this_sac(1),this_sac(2),'k.')
end

xlim([-30,30]);
ylim([-30,30]);