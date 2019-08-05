%% get model comp table
%
% -------------------
% Jeff Mohl
% 8/5/19
% -------------------
%
% Description: given a cell array of model structures, return AIC and BIC
% values across subjects
%

function [AIC,BIC] = get_model_comp_table(m_fits, models, n_params)

AIC = zeros(length(m_fits),length(models));
BIC = zeros(length(m_fits),length(models));

for subj = 1:length(m_fits)
    this_subject = m_fits{subj};
    for model = 1:length(models)
        this_model = models{model};
        model_ind = ismember(vertcat(this_subject.models{:}),this_model,'rows');
        nll = this_subject.nll{model_ind};
        if this_model(3) == 3
            n_obs = sum(this_subject.responses{model_ind}{1,1},'all');
        else
            n_obs = sum(this_subject.responses{model_ind},'all');
        end
        [AIC(subj,model),BIC(subj,model)] = aicbic(-nll,n_params(model),n_obs);
    end
end

% put results in table for easier viewing
% model_comp_table = table({'CI'; 'seg';'int'},[sum(test_nlls.CI);sum(test_nlls.seg);sum(test_nlls.int)],[AIC_CI;AIC_seg;AIC_int], [BIC_CI;BIC_seg;BIC_int],...
%     [AIC_CI;AIC_seg;AIC_int]-best_AIC,[BIC_CI;BIC_seg;BIC_int]-best_BIC, 'VariableNames', {'Model';'nll';'AIC_val'; 'BIC_val';'AIC_diff';'BIC_diff'});
end