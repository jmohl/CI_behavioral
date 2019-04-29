%% Apply post hoc behavioral filtering to monkey data
%
% -------------------
% Jeff Mohl
% 4/25/19
% -------------------
%
% Description: Monkey data has a lot of behavioral quirks that is not
% present in human data. For instance monkeys may lapse mid trial and make
% the second saccade very late, resulting in single saccade trials even
% when the targets are well separated. Also monkeys may make vertical
% saccades which are not appreciably different in the X coordinate, but
% these should be treated as single reports I think.

% NOT FINISHED
function monkey_behavior_filter(data)

%find trials where X coordinate of A and V sacs is less than 3 degrees
%apart.
AV_diff = cell2mat(cellfun(@minus,data.A_endpoints, data.V_endpoints,'UniformOutput',0)); %note that this removes all empty vectors, so to keep the indices need to only apply to 2 sac trials
corrective_ind = abs(AV_diff(:,1)) < 3; %only gives me 81 out of 5502 trials

test = data(data.n_sacs > 1,:);
test_filt = test(corrective_ind,:);

histc(test_filt(:,{'A_tar','V_tar'}))

test = data(data.A_tar == -6 & data.V_tar == -6,:);
test_diff = cell2mat(test.AV_diff);
histogram(test_diff(:,1))
A_mat = cell2mat(test.A_endpoints);

histogram(A_mat(:,1))
figure()
A_alone = data(data.A_tar == -6 & strcmp(data.trial_type,'A'),:).valid_endpoints;
A_alone = cell2mat(A_alone);
histogram(A_alone(:,1));

% find trials where there is only 

end