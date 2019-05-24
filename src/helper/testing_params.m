%for setting params and testing of datalike updates

MAXRNG = 50.5;
fitoptions.binsize = 1; %size of bins used for responses, in degrees.
fitoptions.load_saved_fits = 0; %load saved fits, if they exist
fitoptions.make_plots = 0;
fitoptions.n_iterations = 3; %set option to repeat fminsearch for n times
fitoptions.parameter_names = {'A_sig','V_sig','prior_sig','p_common','lambda_uni','lambda_loc','resp_mu_A','resp_mu_V','resp_sig_A','resp_sig_V'};
fitoptions.UBND = [6 15 40 .9 .25 .25 30 30 40 40]; %upper bounds on theta for grid search
fitoptions.LBND = [.5 1 5 .1 .01 .01 0 0 5 5]; %lower bounds on theta for grid search
fitoptions.grid_fineness = 3; %number of points per parameter in grid search, remember n points in grid = grid_fineness^n_params;; 
fitoptions.fmin_options = optimset('MaxFunEvals',10000,'MaxIter',10000, 'TolFun', 1e-2, 'TolX',1e-2);
fitoptions.eval_range = linspace(-MAXRNG,MAXRNG,MAXRNG*2/fitoptions.binsize + 1);
fitoptions.eval_midpoints = linspace(-MAXRNG+fitoptions.binsize/2,MAXRNG-fitoptions.binsize/2,length(fitoptions.eval_range)-1);
fitoptions.priorloc = [-24 -18 -12 -6 6 12 18 24];

theta = [3 5 15 .5 .2 .01];
prior_type = 2;
unisensory_loc = 1;
eval_midpoints = fitoptions.eval_midpoints;

conds_A = [-24 -6 6 24]';
conds_V = [-24 -18 -12 -6 6 12 18 24]';