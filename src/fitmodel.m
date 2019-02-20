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
% model(1) Model type (1=Bayesian, 2=null, 3=switching)
% model(2) Response type (1=unity judgement, 2 = location)
% model(3) estimation proceedure (1=numerical integration)

%% process data for fitting procedure


function [fit_theta,fit_nll,fit_dist]=fitmodel(conditions,responses,model,fitoptions)
% set relevant variables from options
UBND = fitoptions.UBND;
LBND = fitoptions.LBND;
grid_fineness = fitoptions.grid_fineness;
eval_midpoints = fitoptions.eval_midpoints;
fmin_options = fitoptions.fmin_options;
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
datalike_minsearch = @(theta)datalike(conditions,responses,theta,model,eval_midpoints);

%step 1: evaluate likelihood at all values on grid, pick 5 best points
fprintf('Step 1:grid search\n')
theta_range = zeros(length(UBND),grid_fineness);
for ii = 1:length(UBND)
    theta_range(ii,:) = linspace(UBND(ii),LBND(ii),grid_fineness);
end
[t1,t2,t3,t4,t5] = ndgrid(theta_range(1,:),theta_range(2,:),theta_range(3,:),theta_range(4,:),theta_range(5,:));
grid_like = zeros(size(t1));
for grid_pt = 1:numel(t1)
    this_theta = [t1(grid_pt) t2(grid_pt) t3(grid_pt) t4(grid_pt) t5(grid_pt)];
    grid_like(grid_pt) = datalike_minsearch(this_theta);
end
[~,min_inds] = sort(grid_like(:));
best_thetas =[t1(min_inds(1:5)) t2(min_inds(1:5)) t3(min_inds(1:5)) t4(min_inds(1:5)) t5(min_inds(1:5))];

%step 2: use best grid params as starting point for fminsearch. 
fit_thetas = zeros(size(best_thetas));
fit_nlls = zeros(size(best_thetas,1),1);
for ii = 1:size(best_thetas,1)
    fprintf('Step 2:fmin search, iter:%d\n',ii)
    [fit_thetas(ii,:),fit_nlls(ii),~,~] = fminsearch(datalike_minsearch,best_thetas(ii,:),fmin_options);
end
%this was able to run the optimization procedure in 6.4 sec for the unity
%judgement case. For the localization case it takes much longer.
[fit_nll,best_ind] = min(fit_nlls(:));
fit_theta = fit_thetas(best_ind,:);

[~,fit_dist] = datalike(conditions,responses,fit_theta,model,eval_midpoints);

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








