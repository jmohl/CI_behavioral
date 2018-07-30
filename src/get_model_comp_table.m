%% get model comp table
%
% -------------------
% Jeff Mohl
% 7/30/18
% -------------------
%
% Description: compute AIC and BIC for all models and return a table of
% results comparing them to one another.
%

function model_comp_table = get_model_comp_table(test_nlls, n_params, n_obs)

[AIC_CI,BIC_CI] = aicbic(-test_nlls.CI,n_params.CI,n_obs);
[AIC_seg,BIC_seg] = aicbic(-test_nlls.seg,n_params.seg,n_obs);
[AIC_int,BIC_int] = aicbic(-test_nlls.int,n_params.int,n_obs);

best_BIC = min([BIC_CI,BIC_int,BIC_seg]);
best_AIC = min([AIC_CI,AIC_int,AIC_seg]);

% put resultsi n table for easier viewing
model_comp_table = table({'CI'; 'seg';'int'},[test_nlls.CI;test_nlls.seg;test_nlls.int],[AIC_CI;AIC_seg;AIC_int], [BIC_CI;BIC_seg;BIC_int],...
    [AIC_CI;AIC_seg;AIC_int]-best_AIC,[BIC_CI;BIC_seg;BIC_int]-best_BIC, 'VariableNames', {'Model';'nll';'AIC_val'; 'BIC_val';'AIC_diff';'BIC_diff'});
end