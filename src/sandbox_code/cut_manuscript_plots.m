%% Figure code that didn't make the manuscript but saving just in case
%
% -------------------
% Jeff Mohl
% 8/30/19
% -------------------
%
% Description: These are all figures I originally intended to have in the
% manuscript but decided to cut

%% Condensed localization plot

model = [1 1 3 1];
true_loc = 0; %option to use true target locations or relative locations (from unimodal saccades) for specifying bias.
plot_condensed_loc(models_mj,model,true_loc);
title('Juno')
if savefiles
    saveas(gcf,sprintf('%s\\juno_condensed_bias',figpath),'png');
    saveas(gcf,sprintf('%s\\juno_condensed_bias',figpath),'svg');
end

plot_condensed_loc(models_my,model,true_loc);
title('Yoko')
if savefiles
    saveas(gcf,sprintf('%s\\yoko_condensed_bias',figpath),'png');
    saveas(gcf,sprintf('%s\\yoko_condensed_bias',figpath),'svg');
end

plot_condensed_loc(models_h,model,true_loc);
title('Human')
if savefiles
    saveas(gcf,sprintf('%s\\HU_condensed_bias',figpath),'png');
    saveas(gcf,sprintf('%s\\HU_condensed_bias',figpath),'svg');
end

%% AIC BIC table
%which models to compare
joint_models = {[1 1 3 1];[1 2 3 1];[1 3 3 1]}; %Bayes, model selection, probabilistic fusion (null)
unity_models = {[1 0 1 1];[2 0 1 1]};
loc_models = {[1 1 2 1];[1 2 2 1];[1 3 2 1]};
%comparing models on the localization component only, as the
n_params = [5, 5, 5];

models = unity_models;

[AIC_h,BIC_h] = get_model_comp_table(models_h,models,n_params);
[AIC_mj,BIC_mj] = get_model_comp_table(models_mj,models,n_params);
[AIC_my,BIC_my] = get_model_comp_table(models_my,models,n_params);

%get relative to non-CI model (probabilistic fusion)
AIC_mj_rel = AIC_mj - AIC_mj(:,end);
AIC_my_rel = AIC_my - AIC_my(:,end);
AIC_h_rel = AIC_h - AIC_h(:,end);

%get means + std for all
%rows are [humans;monkey_j;monkey_y]
%columns are [models 1:3]
AIC_mean = [mean(AIC_h_rel,1);mean(AIC_mj_rel,1);mean(AIC_my_rel,1);];
AIC_sem = [std(AIC_h_rel,1)/sqrt(length(models_h));std(AIC_mj_rel,1)/sqrt(length(models_mj));std(AIC_my_rel,1)/sqrt(length(models_my));];
% BIC_mean = [mean(BIC_h_rel,1);mean(BIC_mj_rel,1);mean(BIC_my_rel,1);];
% BIC_sem = [std(BIC_h_rel,1)/sqrt(length(models_h));std(BIC_mj_rel,1)/sqrt(length(models_mj));std(BIC_my_rel,1)/sqrt(length(models_my));];
clear model_names
for ind = 1:length(models)
    model_names{ind} = get_model_names(models{ind});
end
model_comp_AIC = array2table(AIC_mean,'RowNames',{'h','mj','my'},'VariableNames',horzcat(model_names{:}))
% model_comp_BIC = array2table(BIC_mean,'RowNames',{'h','mj','my'},'VariableNames',horzcat(model_names{:}))

%% BIC figure
figure
subj_bars = barh(AIC_mean(:,1:2));
subj_bars(1).FaceColor = [.5 .5 .5];
subj_bars(2).FaceColor = [0 0 0];

yticklabels({'h','mj','my'});
hold on
errorbar(AIC_mean(:,1:2),[.86 1.14;1.86 2.14;2.86 3.14],AIC_sem(:,1:2),'k.','horizontal') %2nd array sets offset so error bars are in the middle of bars

%add individual subjects?
scatter(reshape(AIC_h_rel(:,1:2),1,[]),reshape(repmat([.86 1.14],7,1),1,[]),[],[1:7,1:7],'filled');
scatter(reshape(AIC_mj_rel(:,1:2),1,[]),reshape(repmat([1.86 2.14],10,1),1,[]),[],[1:10,1:10],'filled');
scatter(reshape(AIC_my_rel(:,1:2),1,[]),reshape(repmat([2.86 3.14],10,1),1,[]),[],[1:10,1:10],'filled');

legend(subj_bars,horzcat(model_names{1:2}),'Interpreter','none')
ylabel('Subjects')
xlabel('Delta AIC - rel non-CI')
if savefiles
    saveas(gcf,sprintf('%s\\deltaAIC_loc',figpath),'png');
    saveas(gcf,sprintf('%s\\deltaAIC_loc',figpath),'svg');
end
