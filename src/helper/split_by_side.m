%splitting data into left and right subsets

subject = 'Yoko';

load(sprintf('data\\%s_combined.mat',subject));

left_data = tidy_data(tidy_data.A_tar < 0 | tidy_data.V_tar < 0,:); 
%remove conditions that also include a target on the right
left_data = left_data(~(left_data.A_tar > 0  | left_data.V_tar > 0),:); 

right_data = tidy_data(tidy_data.A_tar > 0 | tidy_data.V_tar > 0,:); 
%remove conditions that also include a target on the right
right_data = right_data(~(right_data.A_tar < 0  | right_data.V_tar < 0),:); 

%want the variable to be named tidy_data so I don't mess up my code
tidy_data = left_data;
save(sprintf('data\\%s_left_combined',subject),'tidy_data');
tidy_data = right_data;
save(sprintf('data\\%s_right_combined',subject),'tidy_data');
