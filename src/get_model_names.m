%% get string model names, from model numbers
%
% -------------------
% Jeff Mohl
% 3/8/19
% -------------------
%
% Description: convert from numerical model description used in fitting procedure to strings.
function model_names = get_model_names(models)
model_names = cell(1,size(models,2));
for i = 1:size(models,2)
    model_str = [];
    switch models{1,i}(1)
        case 2
            model_str = 'B';
        case 3
            model_str = 'F';
        case 4
            model_str = 'MS';
    end
    switch models{1,i}(2)
        case 1
            model_str = strcat(model_str,'_U'); %unity judgement
        case 2
            model_str = strcat(model_str,'_L'); %location judgement
    end
    model_names{i} = model_str;
end
    