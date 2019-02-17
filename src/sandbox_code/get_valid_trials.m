%% Trial filtration
%
% -------------------
% Jeff Mohl
% 2/17/19
% -------------------
%
% Description: Filtering trials more aggressively for behavioral analysis.
% Trials will only be included if they meet the following criteria.
% 1. Go cue was reached (trial was initiated)
% 2. duration of trial was at least 700 ms after go cue (removes trials
% where trial ended because of saccade inaccuracy before all reports could be
% made)
% 3. ...

% random helper code for building function

addpath('C:\Users\jtm47\Documents\Projects\AVD_data_exploration\src')

data = load('Juno_AVD2_2017_07_31_2_tidy.mat');
data = data.tidy_data;

function valid_data = get_valid_trials(data)

valid_data = data(~isnan(data.go_time)&(data.end_time - data.go_time > 700),:);


% A_tar = 6;
% V_tar = 12;
% this_data = data(data.A_tar == A_tar & data.V_tar == V_tar,:);
% this_valid_data = valid_data(valid_data.A_tar == A_tar & valid_data.V_tar == V_tar,:);
% plot_eye_traces(this_data);
% plot_eye_traces(this_valid_data);
end