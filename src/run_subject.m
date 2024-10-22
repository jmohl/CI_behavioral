%% script for running a single daaset through the analysis
%
% -------------------
% Jeff Mohl
% 2/15/19
% -------------------
%
% Description: this script will allow running of each step of the
% behavioral CI project, from loading data, fitting model, and generating
% figures, for a single subject. It is meant to be run out of
% master_script and so at least the hard coded section of that script is
% assumed to have been run. UPDATED to run on new modelfit procedure.

%note: H01 has some problems and it is best to exclude for now, since I
%will probably have to rerun that subject in any case

function run_subject(subject,raw_data,model_list)
global fitoptions

%% Clean data

%get valid data only
valid_data = raw_data(logical(raw_data.valid_tr),:);

%omit trials without valid endpoints, very rarely removes data
data = valid_data(~cellfun('isempty',valid_data.valid_endpoints),:);


%% split data into train and test sets
if fitoptions.cross_validate
[AV_train,AV_test] = get_train_test(data,fitoptions.kfolds);
end
%% CI model (new)
clear m;
%load saved modelfits, if exist
if exist(sprintf('results\\modelfits\\%s_m.mat',subject),'file') 
    m=load(sprintf('results\\modelfits\\%s_m.mat',subject));
    m=m.m;
else % if there is no modelstruct at all, initialize model struct
    m.fitoptions = fitoptions;
    m.subject = subject;
    m.models={};
    m.nll={};
    m.thetas = {};
    m.fit_dist={};
    m.conditions={};
    m.responses={};
end

for mi = 1:size(model_list,1)
    if exist(sprintf('results\\modelfits\\%s_m.mat',subject),'file')
        model_exists = ismember(model_list{mi},vertcat(m.models{:}),'rows');
    else
        model_exists = 0;
    end
    if ~isempty(m.models)&& model_exists && fitoptions.load_saved_fits
        fprintf('skipping model [%d %d %d %d] because already saved \n',model_list{mi})
        continue
    end
    if fitoptions.cross_validate % if using cross validation
        fit_theta = cell(fitoptions.kfolds,1);
        fit_dist = cell(fitoptions.kfolds,1);
        fit_nll = cell(fitoptions.kfolds,1);
        model = model_list{mi};
        
        for ki = 1:fitoptions.kfolds
            train_data = AV_train{ki};
            test_data = AV_test{ki};
            fprintf('Fitting Subject: %s, Model: %d %d %d %d, k-fold:%d\n',subject,model, ki)
            [conditions,responses] = get_prepro_data(train_data,model);
            [fit_theta{ki},~,fit_dist{ki}]=fitmodel(conditions,responses,model);
            [conditions,responses] = get_prepro_data(test_data,model);
            fit_nll{ki} = datalike(conditions,responses,fit_theta{ki},model,fitoptions.eval_midpoints);
        end
         %for conditions and responses (used for plotting), return whole dataset
            [conditions,responses] = get_prepro_data(data,model);
    else
        model = model_list{mi}; 
        fprintf('Fitting Subject: %s, Model: %d %d %d %d\n',subject,model)
        [conditions,responses] = get_prepro_data(data,model);
        if fitoptions.use_uni_means % testing JM
            mean_locs = get_unimodal_means(conditions,data,model);
            [fit_theta,fit_nll,fit_dist]=fitmodel(mean_locs,responses,model);  
        else
            [fit_theta,fit_nll,fit_dist]=fitmodel(conditions,responses,model);
        end

    end
    %if the model structure already has the given model, then overwrite,
    %otherwise add new row
    if model_exists 
        model_ind = ismember(vertcat(m.models{:}),model_list{mi},'rows');
        m.models{model_ind} = model;
        m.thetas{model_ind} = fit_theta;
        m.nll{model_ind} = fit_nll;
        m.fit_dist{model_ind} = fit_dist;
        m.conditions{model_ind} = conditions;
        m.responses{model_ind} = responses;
    else
        m.models{end+1} = model;
        m.thetas{end+1} = fit_theta;
        m.nll{end+1} = fit_nll;
        m.fit_dist{end+1} = fit_dist;
        m.conditions{end+1} = conditions;
        m.responses{end+1} = responses;
    end
    %save out model struct
    save(sprintf('results\\modelfits\\%s_m',subject),'m');
end
clear m rawdata; %not convinced this does anything

end