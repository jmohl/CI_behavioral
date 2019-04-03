%% get string model names, from model numbers
%
% -------------------
% Jeff Mohl
% 3/8/19
% -------------------
%
% Description: convert from numerical model description used in fitting procedure to strings.
function model_names = get_model_names(models)
model_names = cell(1,size(models,1));
for i = 1:size(models,1)
    model_str = [];
    switch models(i,1) %unity judgement type
        case 1
            model_str = 'B'; 
        case 2
            model_str = 'PF';
    end
    switch models(i,3) %location judgement type
        case 1
            model_str = strcat(model_str,'_B');
        case 2
            model_str = strcat(model_str,'_MS');
        case 3
            model_str = strcat(model_str,'_PF');
    end
    switch models(i,2)  %fit type
        case 1
            model_str = strcat(model_str,'_U'); %unity judgement
        case 2
            model_str = strcat(model_str,'_L'); %location judgement
        case 3
            model_str = strcat(model_str,'_J'); %joint fit
    end
    model_names{i} = model_str;
end
    