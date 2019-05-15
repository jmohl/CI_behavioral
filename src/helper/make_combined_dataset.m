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
%dates used     
%     {'Juno_AVD2_2017_07_20_2'}
%     {'Juno_AVD2_2017_07_26'  }
%     {'Juno_AVD2_2017_08_02'  }
%     {'Juno_AVD2_2017_09_06'  }
%     {'Juno_AVD2_2017_12_13'  }
%     {'Juno_AVD2_2018_01_23'  }
%     {'Juno_AVD2_2018_02_27'  }
%     {'Juno_AVD2_2018_04_05'  }
%     {'Juno_AVD2_2018_04_09'  }
%     {'Juno_AVD2_2018_04_10'  }
% 6476 valid trials, 11128 total trials, for 10 pooled days.
save('data\Juno_combined','tidy_data');

subject = 'Yoko';
n_pooled_days = 10;
seed = 'default';

tidy_data = load_pool_data(n_pooled_days,subject,seed);
%     {'Yoko_AVD2_2018_11_29'}
%     {'Yoko_AVD2_2019_01_09'}
%     {'Yoko_AVD2_2019_01_10'}
%     {'Yoko_AVD2_2019_01_30'}
%     {'Yoko_AVD2_2019_03_07'}
%     {'Yoko_AVD2_2019_03_12'}
%     {'Yoko_AVD2_2019_03_19'}
%     {'Yoko_AVD2_2019_04_03'}
%     {'Yoko_AVD2_2019_04_17'}
%     {'Yoko_AVD2_2019_04_24'}
% 5830 valid trials, 8673 total trials, for 10 pooled daya
save('data\Yoko_combined','tidy_data');
