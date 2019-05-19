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
% set relevant variables from options
UBND = fitoptions.UBND;
LBND = fitoptions.LBND;
debug = 0;

%% start of actual fitting procedure
%adapting from Acerbi, this fitting procedure will progress in 2 steps to
%start. First, a grid spanning reasonable paramter space will be evaluated,
%and the best ## points will be chosen and used as initial theta values for
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

n_iter = fitoptions.n_iterations;
%step 1: evaluate likelihood at all values on grid, pick 5 best points
fprintf('Step 1:grid search\n')
tic
best_thetas = grid_search(datalike_minsearch);
best_thetas = best_thetas(1:n_iter,:); %pick only the desired number to iterate on
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

if debug 
    %plot some things for comparing with behavior
    if model(2) == 1
        plot_psingle(responses,conditions,fit_dist);
    end
    if model(2) == 2
        for ic = 1:length(conditions)
            %plotting the real saccade distributions and those predicted by
            %the model for every condition
            plot_modelhist(responses(ic,:),fit_dist(ic,:),eval_midpoints)
            title(sprintf('%d A %d V',conditions(ic,:)))
            set(gcf,'Position',[100,60,1049,895])
        end
    end
        
end

end








