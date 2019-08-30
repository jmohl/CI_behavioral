%% Fit model to data
%
% -------------------
% Jeff Mohl
% 2/1/19
% -------------------
%
% Description: this is code for performing the model fitting steps for the
% given model and given subject. Will return the best fit parameters, the
% nll for the model fit... maybe more in the future
%
% Inputs:
% data(1) trial number
% data(2) A tar
% data(3) V tar
% data(4) response (number of saccades for unity judgement case)
% theta - contains model fit parameters  
% model(1) Model type (1=Bayesian, 2=updated bayesian, 3=naive)
% model(2) Response type (1=unity judgement, 2 = location, 3= joint fit)
% model(3) estimation proceedure (1=numerical integration)

%% process data for fitting procedure


function [fit_theta,fit_nll,fit_dist]=fitmodel(conditions,responses,model)

global fitoptions
%% start of actual fitting procedure
%adapting from Acerbi, this fitting procedure will progress in 2 steps to
%start. First, a grid spanning reasonable paramter space will be evaluated,
%and the best n_iterations points will be chosen and used as initial theta values for
%fminsearch. In the future might switch from fminsearch to bads if that
%works better.

%todos: validation

%2/1/19 note: used timeit to test this likelihood function as 0.0091 sec
%for human subject H08 unity judgement. Unfortunately for the localization
%judgement case it takes 2.28 seconds (with eval_n = 250, with eval_n = 100
%significantly faster at 0.147 sec), which is pretty slow considering the
%number of times this function needs to be called. Also the likelihood does
%not converge over 500 iterations, which takes 2 minutes.
%
% f=@() datalike_minsearch(theta);
% timeit(f)
datalike_minsearch = @(theta)datalike(conditions,responses,theta,model,fitoptions.eval_midpoints);

%step 1: evaluate likelihood at all values on grid, pick n_iteration best points
fprintf('Step 1:grid search\n')
[UBND,LBND] = get_ini_params(model); %inital points are hardcoded below
tic
best_thetas = grid_search(datalike_minsearch,UBND,LBND,fitoptions.grid_fineness);
best_thetas = best_thetas(1:fitoptions.n_iterations,:); %pick only the desired number to iterate on
toc
%step 2: use best grid params as starting point for fminsearch. 
fit_thetas = zeros(size(best_thetas));
fit_nlls = zeros(size(best_thetas,1),1);
for ii = 1:size(best_thetas,1)
    fprintf('Step 2:fmin search, iter:%d\n',ii)
    tic
    [fit_thetas(ii,:),fit_nlls(ii),~,~] = fminsearch(datalike_minsearch,best_thetas(ii,:),fitoptions.fmin_options);
    toc
end
%this was able to run the optimization procedure in 6.4 sec for the unity
%judgement case. For the localization case it takes much longer.
[fit_nll,best_ind] = min(fit_nlls(:));
fit_theta = fit_thetas(best_ind,:);

[~,fit_dist] = datalike(conditions,responses,fit_theta,model,fitoptions.eval_midpoints);

end

function [UBND,LBND] = get_ini_params(model)

% model values
% model(1) = CI type: Bayes (1) or probabilistic fusion (2)
% model(2) = stimulus fusion: Bayes reweight(1), model selection(2), probabilistic fusion (3)
% model(3) = task/fit type: unity judgement (1), localization (2), joint fit (3), unisensory localization (4)
% model(4) = prior type: naive normal (1), discrete empirical (2), normal mixture empirical (3)

% potential theta values
% Included in all models:
% theta(1) = Visual sigma
% theta(2) = auditory sigma
% theta(3) = prior sigma, variance of prior (not used in discrete prior case)
% theta(4) = p_common, prior probability of common cause
% theta(5) = lapse probability on unity judgement
% Optional thetas
% theta(6) = prior_mu, for when using a mixture of normals prior
% theta(7) prior_mu2
% theta(8) prior_sig2

UBND = [6 15 40 .9 .25];
LBND = [.5 1 5 .1 0.01];


if model(4) == 3 
    UBND = [UBND,15,15];%,15,25];
    LBND = [LBND,-15,5];%,-15,5];
end

end






