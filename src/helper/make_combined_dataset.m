%% pooling data from 10 random days, to save as a single dataset
%
% -------------------
% Jeff Mohl
% 3/4/19
% -------------------
%
% Description: was previously pooling different data for each model fit
% run, but that led to inconsistent results and made model comparison
% impossible. Instead here I am pooling days once, and saving that as a new
% data file for use in the model.

% note that the folder needs to be set to the CI_behavioral project for
% this to run properly.

%pool for juno
subject = 'Juno';
n_pooled_days = 10;
seed = 'default';

tidy_data = load_pool_data(n_pooled_days,subject,seed);
%dates used {'Juno_AVD2_2017_07_11';'Juno_AVD2_2017_07_13';'Juno_AVD2_2017_07_19';
%   'Juno_AVD2_2017_08_30';'Juno_AVD2_2017_12_07';'Juno_AVD2_2018_01_05';
%   'Juno_AVD2_2018_02_21_2';'Juno_AVD2_2018_03_19_2';'Juno_AVD2_2018_04_09';
%   'Juno_AVD2_2018_04_10'}
% 6749 valid trials, 11443 total trials
save('data\Juno_combined','tidy_data');

subject = 'Yoko';
n_pooled_days = 10;
seed = 'default';

tidy_data = load_pool_data(n_pooled_days,subject,seed);
% {'Yoko_AVD2_2018_11_16';'Yoko_AVD2_2018_11_20';'Yoko_AVD2_2018_11_27';
% 'Yoko_AVD2_2018_11_29';'Yoko_AVD2_2018_12_13';'Yoko_AVD2_2019_01_09';
% 'Yoko_AVD2_2019_01_10';'Yoko_AVD2_2019_01_11';'Yoko_AVD2_2019_01_18';
% 'Yoko_AVD2_2019_01_31'}
% 5582 valid trials, 7749 total trials
save('data\Yoko_combined','tidy_data');
