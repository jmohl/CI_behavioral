%% removing spike data files
%
% -------------------
% Jeff Mohl
% 10/22/19
% -------------------
%
% Description: removing spike data from data files in this project
%

 local_directory = 'C:\Users\jtm47\Documents\Projects\CI_behavioral\data';
 cd(local_directory)

juno_data = dir('*Juno*');
yoko_data = dir('*Yoko*');

for y_ind = 1:length(yoko_data)
    this_name = yoko_data(y_ind).name;
    this_data = load(this_name);
    tidy_data = this_data.tidy_data;
    if ismember('spkdata', tidy_data.Properties.VariableNames)
        tidy_data = removevars(tidy_data,{'spkdata'});
    end
    save(this_name,'tidy_data');
end


for j_ind = 1:length(juno_data)
    this_name = juno_data(j_ind).name;
    this_data = load(this_name);
    tidy_data = this_data.tidy_data;
    if ismember('spkdata', tidy_data.Properties.VariableNames)
        tidy_data = removevars(tidy_data,{'spkdata'});
    end
    save(this_name,'tidy_data');
end