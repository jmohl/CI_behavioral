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
if exist(sprintf('results\\modelfits\\%s_m.mat',subject),'file') &&  fitoptions.load_saved_fits
    m=load(sprintf('results\\modelfits\\%s_m.mat',subject));
    m=m.m;
else %initialize model struct
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
    if ~isempty(m.models)&& ismember(model_list{mi},vertcat(m.models{:}),'rows')
        fprintf('skipping model [%d %d %d %d] because already saved \n',model_list{mi})
        continue
    end
    if fitoptions.cross_validate % if using cross validation
        fit_theta = cell(fitoptions.kfolds,1);
        fit_dist = cell(fitoptions.kfolds,1);
        test_nll = cell(fitoptions.kfolds,1);
        model = model_list{mi};
        
        for ki = 1:fitoptions.kfolds
            train_data = AV_train{ki};
            test_data = AV_test{ki};
            fprintf('Fitting Subject: %s, Model: %d %d %d %d, k-fold:%d\n',subject,model, ki)
            [conditions,responses] = get_prepro_data(train_data,model);
            [fit_theta{ki},~,fit_dist{ki}]=fitmodel(conditions,responses,model);
            [conditions,responses] = get_prepro_data(test_data,model);
            test_nll{ki} = datalike(conditions,responses,fit_theta{ki},model,fitoptions.eval_midpoints);
        end
            m.models{end+1} = model;
            m.thetas{end+1} = fit_theta;
            m.nll{end+1} = test_nll;
            m.fit_dist{end+1} = fit_dist;
            [m.conditions{end+1},m.responses{end+1}] = get_prepro_data(data,model); %for conditions and responses (used for plotting), return whole dataset
    else
        model = model_list{mi}; 
        fprintf('Fitting Subject: %s, Model: %d %d %d %d\n',subject,model)
        m.models{end+1} = model;
        [conditions,responses] = get_prepro_data(data,model);
        [fit_theta,fit_nll,fit_dist]=fitmodel(conditions,responses,model);
        m.thetas{end+1} = fit_theta;
        m.nll{end+1} = fit_nll;
        m.fit_dist{end+1} = fit_dist;
        m.conditions{end+1} = conditions;
        m.responses{end+1} = responses;
    end
    
    %save out model struct
    save(sprintf('results\\modelfits\\%s_m',subject),'m');
end
clear m rawdata; %not convinced this does anything but am getting a lot of crashing

%% making plots for single subject
%some evaluation plots, note that this replots ALL models in the current
%configuration. will probably change that in the future but is fine to just
%set make_plots = 0 for now
if fitoptions.make_plots
    set(0,'DefaultFigureVisible','off');
    for mi = 1:size(model_list,1)
        model = m.models{mi};
        
        if ismember(model(2), [1 3])
            try
                mkdir(sprintf('results\\p_single\\model%d%d%d',model))
            end
            if model(2) == 3
                fit_dist = m.fit_dist{mi}{1};
                responses = m.responses{mi}{1};
            else 
                fit_dist = m.fit_dist{mi};
                responses = m.responses{mi};
            end
            plot_psingle(responses,m.conditions{mi},fit_dist);
            saveas(gcf,sprintf('results\\p_single\\model%d%d%d\\psing_%s',model,m.subject),'png');
        end
        
        if ismember(model(2), [2 3])
            try
                mkdir(sprintf('results\\localization\\model%d%d%d\\%s',model,m.subject))
            end
            for ic = 1:length(m.conditions{mi})
                %plotting the real saccade distributions and those predicted by
                %the model for every condition ic = 5;if model(2) == 3
                if model(2) == 3
                    fit_dist = m.fit_dist{mi}{2};
                    responses = m.responses{mi}{2};
                else
                    fit_dist = m.fit_dist{mi};
                    responses = m.responses{mi};
                end
                plot_modelhist(responses(ic,:,:),fit_dist(ic,:,:),m.fitoptions.eval_midpoints)
                title(sprintf('%d A %d V',m.conditions{mi}(ic,:)))
                set(gcf,'Position',[100,60,1049,895])
                saveas(gcf,sprintf('results\\localization\\model%d%d%d\\%s\\%dA%dV',model,m.subject,m.conditions{mi}(ic,:)),'png');
            end
        end
    end
    close all;
    set(0,'DefaultFigureVisible','on');
end
end