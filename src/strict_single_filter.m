%% working on more stringent behavioral filtering criteria
%
% -------------------
% Jeff Mohl
% 7/3/19
% -------------------
%
% Description: some subjects, especially yoko, like to make single saccades
% that are to one of two well separated targets (essentially ignoring the
% task). Here I'm going to try and filter out this kind of trial
% specifically. Also making plots to compare the filtered and unfiltered
% data.

%endpoints are already separated out. Specifically I want single endpoints
%from conditions with at least 12 degrees of target sep.
function tidy_data = strict_single_filter(tidy_data)
rel_trs = abs(tidy_data.A_tar - tidy_data.V_tar) >= 12 & tidy_data.n_sacs == 1;

just_horz = cell2mat(tidy_data(rel_trs,:).valid_endpoints);
just_horz = just_horz(:,1);

V_error = just_horz - tidy_data(rel_trs,:).V_tar;
A_error = just_horz - tidy_data(rel_trs,:).A_tar;

exclusions = abs(V_error) > 10 | abs(A_error) > 20;
tidy_data(rel_trs,:).valid_tr(exclusions) = 0;

% histogram(V_error)
%from this histogram you can clearly see all the outlier trials. That is
%pretty much exactly what I'm talking about. There are 125 of these trials
%with error > 10, out of 1122 trials that I am considering.
% histogram(A_error)
%this histogram however is really interesting. It's super bio-modally
%distributed centered around ~15 deg error. This makes it seem like she is
%much more likely to make just a visual saccade in conditions where the two
%targets are close together (12 degrees, maybe even 18 degrees) This makes
%me think that the filter for this condition should be > 20 degrees, which
%is still 86 trials.

%if I want to be more flexible I can calculate the normal error
%distribution including the A and V sacs, and exlude any outlier trials
%that way. For now will just look at excluding trials using this fixed
%criteria. 
end