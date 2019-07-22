%% Collapsing across target disparities
%
% -------------------
% Jeff Mohl
% 7/19/19
% -------------------
%
% Description: collapsing localization information across multiple target
% pairs, to make a single plot that captures all information and can be
% compared between species
%

function plot_condensed_loc(m,model,true_loc)
%starting with a single subjects responses
if length(m) > 1
   % pool data across model structures if provided with multiple 
   % otherwise, get data from single model structure
   % responses and fit distributions are 3d, so add another dimension to go
   % across subjects then take mean across this dimension.
   for subj_ind = 1:length(m)
       this_m = m{subj_ind};
       model_ind = ismember(vertcat(this_m.models{:}),model,'rows');
       conditions = this_m.conditions{model_ind};
       xrange = this_m.fitoptions.eval_midpoints;
       try
       if model(3) == 3 %joint fit models have cell array rather than vectors for these
           responses(:,:,:,subj_ind) = this_m.responses{model_ind}{2};
           fit_dist(:,:,:,subj_ind) = this_m.fit_dist{model_ind}{2};
       else
           responses(:,:,:,subj_ind) = this_m.responses{model_ind};
           fit_dist(:,:,:,subj_ind) = this_m.fit_dist{model_ind};
       end
       %get unimodal responses
       uni_ind = ismember(vertcat(this_m.models{:}),[0 0 4 1],'rows');
       uni_cond = this_m.conditions{uni_ind};
       A_resp(:,:,:, subj_ind) = this_m.responses{uni_ind}{1};
       V_resp(:,:,:, subj_ind) = this_m.responses{uni_ind}{2};
       end
   end
   %get mean model and response fits
    responses = mean(responses,4);
    fit_dist = mean(fit_dist,4);
    A_resp = mean(A_resp,4);
    V_resp = mean(V_resp,4);
else
    %otherwise, get data from single model structure
    model_ind = ismember(vertcat(m.models{:}),model,'rows');
    conditions = m.conditions{model_ind};
    xrange = m.fitoptions.eval_midpoints;
    if model(3) == 3 %joint fit models have cell array rather than vectors for these
        responses = m.responses{model_ind}{2};
        fit_dist = m.fit_dist{model_ind}{2};
    else
        responses = m.responses{model_ind};
        fit_dist = m.fit_dist{model_ind};
    end
    %get unimodal responses
    uni_ind = ismember(vertcat(m.models{:}),[0 0 4 1],'rows');
    uni_cond = m.conditions{uni_ind};
    A_resp = m.responses{uni_ind}{1};
    V_resp = m.responses{uni_ind}{2};
end

%split localizations into auditory, visual, and combined pdfs
diag_array = diag(ones(length(xrange),1));
diag_array = logical(diag_array);

cA = squeeze(sum(responses(:,:,:),2));%marginalize across V sac
cV = squeeze(sum(responses(:,:,:),3));%marginalize across A sac
cAV = squeeze(responses(:,diag_array));%get only along diagonal
cA_uni= squeeze(A_resp(:,:,:,:));
cV_uni = squeeze(V_resp(:,:,:,:));
%separate single and dual saccades
cA = cA - cAV;
cV = cV - cAV;

%interpolation code, slightly more accurate than using the max
respV_c2 = sum(cV.*xrange,2)./sum(cV,2);
respA_c2 = sum(cA.*xrange,2)./sum(cA,2);
respAV_c1 = sum(cAV.*xrange,2)./sum(cAV,2);
respA_uni = sum(cA_uni.*xrange,2)./sum(cA_uni,2);
respV_uni  = sum(cV_uni.*xrange,2)./sum(cV_uni,2);

A_tar = conditions(:,1);
V_tar = conditions(:,2);
A_cond = uni_cond{1};
V_cond = uni_cond{2};
%have to repeat the mean values for all the different combinations
for tar_ind = 1:length(A_tar)
    respA(tar_ind,1) =  respA_uni(A_cond == A_tar(tar_ind));
end
for tar_ind = 1:length(V_tar)
    respV(tar_ind,1) =  respV_uni(V_cond == V_tar(tar_ind));
end

if true_loc
    AV_dif = A_tar - V_tar;
    A_bias = respA_c2 - A_tar;
    V_bias = respV_c2 - V_tar ;
    AV_bias_A = respAV_c1 - A_tar;
else
    AV_dif = A_tar - V_tar;
    A_bias = respA_c2 - respA;
    V_bias = respV_c2 - respV ;
    AV_bias_A = respAV_c1 - respA;
end

%group by target separation
[gav,gav_labels] = findgroups(AV_dif);
mean_A_bias = splitapply(@mean, A_bias(:,1),gav);
mean_V_bias = splitapply(@mean, V_bias(:,1),gav);
mean_AV_bias = splitapply(@mean, AV_bias_A(:,1),gav);

%plotting relative to true target locations
figure
hold on
%plotting individual points

plot(gav_labels,mean_A_bias,'r-')
plot(gav_labels,mean_V_bias,'b-')
plot(gav_labels,mean_AV_bias,'k-')
plot([25 -25], [-25 25], '--','Color',[.75 .75 .75 .2])
xlabel({'AV (A - V)';'<<<< A target to left of V | A target to right of V >>>>'})
ylabel({'Bias from target location';'(resp - target)'})
legend('A Saccade','V saccade','Single Saccade relative to A tar','ref line complete visual dominance','location', 'best')
plot(AV_dif,A_bias,'.r')
plot(AV_dif,V_bias,'.b')
plot(AV_dif,AV_bias_A,'.k')
title('Response bias vs target separation')

end
