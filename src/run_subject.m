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

%todo: put this in function form to reduce workspace clutter

%% Clean data
%get valid data that reached go time and had an appropriate trial duration
valid_data = raw_data(~isnan(raw_data.go_time)&(raw_data.end_time - raw_data.go_time > 700),:);

% add valid endpoints field that extracts endpoints that are defines as
% "responses" as well as I can define that.
[valid_endpoints,A_endpoints,V_endpoints] = get_response_endpoints(valid_data,0,100);
valid_data.valid_endpoints = valid_endpoints;
valid_data.A_endpoints = A_endpoints;
valid_data.V_endpoints = V_endpoints;

%replace n_sacs field with one calculated using valid response endpoints;
valid_data.n_sacs = cellfun(@(x) size(x,1),valid_data.valid_endpoints);

%omit trials without valid endpoints, very rarely removes data
data = valid_data(~cellfun('isempty',valid_data.valid_endpoints),:);

%% TODO correct eye tracker calibration
if fitoptions.correct_bias 
    [data] = get_bias_corrected_data(data); 
end

%% split data into train and test sets TODO
%[AV_train,AV_test] = get_train_test(data,CI_opts.k_folds,fixed_params.AV_pairs);

%% CI model (new)
clear m;
%load saved modelfits, if exist
if exist(sprintf('results\\modelfits\\%s_m.mat',subject),'file') &&  fitoptions.load_saved_fits
    m=load(sprintf('results\\modelfits\\%s_m.mat',subject));
    m=m.m;
else %initialize model struct
    m.fitoptions = fitoptions;
    m.subject = subject;
end

for mi = 1:size(model_list,1)
    if isfield(m,'models')&& ismember(model_list{mi},vertcat(m.models{:}),'rows')
        fprintf('skipping model [%d %d %d] because already saved \n',model_list{mi})
        continue
    end
    model = model_list{mi}; 
    fprintf('Fitting Subject: %s, Model: %d %d %d\n',subject,model)
    m.models{mi} = model;
    [conditions,responses] = get_prepro_data(data,model,fitoptions);
    [fit_theta,fit_nll,fit_dist]=fitmodel(conditions,responses,model,fitoptions);
    m.thetas{mi} = fit_theta;
    m.nll{mi} = fit_nll;
    m.fit_dist{mi} = fit_dist;
    m.conditions{mi} = conditions;
    m.responses{mi} = responses;
    
    %save out model struct
    save(sprintf('results\\modelfits\\%s_m',subject),'m');
end

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
                mkdir(sprintf('results\\p_single\\model%d%d%d\\%s',model,m.subject))
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
%% START OF OLD CODE
%% fit CI model
% 
% % run fminsearch to fit parameters
% for this_fold = 1:CI_opts.k_folds
% CI_minsearch = @(free_param_vector)get_nll_CI(AV_train(this_fold,:),fixed_params,free_param_vector);
% [fit_results.CI_fit(this_fold,:),fit_results.CI_nll,~,~] = fminsearch(CI_minsearch,param_vector,fmin_options);
% 
% %get nll on test dataset, using fit params
% test_nlls.CI(this_fold) = get_nll_CI(AV_test(this_fold,:),fixed_params,fit_results.CI_fit(this_fold,:));
% end
% 
% %saving out average fit params for plotting later
% fit_params(2,:) = mean(fit_results.CI_fit(:,:));
% 
% %% Alternative models
% %alternative models for comparing to full CI model
% 
% % fully segregated model, 4 params
% param_vector = param_vector(1:4);
% for this_fold = 1:CI_opts.k_folds
% seg_minsearch = @(param_vector)get_nll_seg(AV_train(this_fold,:),fixed_params,param_vector);
% [fit_results.seg_fit,fit_results.seg_nll,~,~] = fminsearch(seg_minsearch,param_vector,fmin_options);
% 
% %get nll on test dataset, using fit params
% test_nlls.seg(this_fold) = get_nll_seg(AV_test(this_fold,:),fixed_params,fit_results.seg_fit);
% 
% % fully integrated model, 4 params
% param_vector = param_vector(1:4);
% int_minsearch = @(param_vector)get_nll_int(AV_train(this_fold,:),fixed_params,param_vector);
% [fit_results.int_fit,fit_results.int_nll,~,~] = fminsearch(int_minsearch,param_vector,fmin_options);
% 
% %get nll on test dataset, using fit params
% test_nlls.int(this_fold) = get_nll_int(AV_test(this_fold,:),fixed_params,fit_results.int_fit);
% 
% end
% fit_params(3,1:4) = mean(fit_results.seg_fit(:,:));
% fit_params(4,1:4) = mean(fit_results.int_fit(:,:));
% 
% 
% %% probability matching instead of posterior reweighted
% % todo
% 
% %% CI model with target locs and sigmas determined by unimodal fits, 2 params
% % get single modality estimates (todo)
% % [fixed_params.V_mu,fixed_params.V_sig] = get_unimodal_est(data.V);
% % [fixed_params.A_mu,fixed_params.A_sig] = get_unimodal_est(data.A);
% % free_params_vector_fAV = free_params_vector(4:end);
% % 
% % CI_fAV_minsearch = @(free_params_vector_fAV)get_nll_CI_fixedAV(AV_tar_pairs,AV_endpoint_train,fixed_params,free_params_vector_fAV);
% % [fit_params_vector_fAV,min_nll_CI_fAV,~,~] = fminsearch(CI_fAV_minsearch,free_params_vector_fAV,options);
%  
% %notes: model does not converge, ends up with extremely large prior sig and
% %extremely low p_common indicating that these variables are adding very
% %little. expectation is that the lack of accuracy is really hurting this
% %model and the parameters are trying to account for that somehow. 
% 
% %% convert fit results to table
% fit_params = array2table(fit_params);
% fit_params.Properties.VariableNames = fieldnames(free_params);
% fit_params.Properties.RowNames = {'initial', 'CI', 'Seg','Int'};
% try 
%     mkdir('results\fit_params')
% end
% save(sprintf('results\\fit_params\\fit_params_%s.mat',subject),'fit_params');
% 
% 
% %% model comparison
% n_params.CI = 5;
% n_params.seg = 4;
% n_params.int = 4;
% 
% % put results in table for easier viewing
% subject
% model_comp_table = get_model_comp_table(test_nlls,n_params,length(vertcat(AV_test{:}))) %n obs = total number of valid saccades in dual trials
% save(sprintf('results\\model_comp\\mc_%s.mat',subject),'model_comp_table');
% writetable(model_comp_table,sprintf('results\\model_comp\\mc_%s.csv',subject));
% 
