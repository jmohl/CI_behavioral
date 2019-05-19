%% run an unbiased grid search on initial parameters
%
% -------------------
% Jeff Mohl
% 5/19/19
% -------------------
%
% Description: takes in the initial ranges for theta, set in master_script,
% and uses them to run an unbiased grid search over the entire parameter
% space. This can get VERY computationally taxing, depending on the number
% of parameters and the fineness of the grid (number of eval points = grid_fineness^n_params)

%inputs: datalike_minsearch - function for evaluating likelihood under a
%given theta, set in fitmodel code.

%outputs: best_thetas - each row is a set of values for theta, with the
%best row followed by the second best etc for all combinations of initial
%parameters.

function best_thetas = grid_search(datalike_minsearch)
global fitoptions

%make grid
for ii = 1:length(fitoptions.UBND)
    theta_range{ii} = linspace(fitoptions.UBND(ii),fitoptions.LBND(ii),fitoptions.grid_fineness);
end
[theta_pts{1:length(fitoptions.UBND)}] = ndgrid(theta_range{:});

%get likelihood at each grid point
 grid_like = arrayfun(@(varargin) datalike_minsearch(cell2mat(varargin)), theta_pts{:});
%get the best parameter values, in order
[~,min_inds] = sort(grid_like(:));
find_best_values = @(ind,min_inds) theta_pts{ind}(min_inds);
for this_param = 1:length(fitoptions.UBND)
    best_thetas(:,this_param) = find_best_values(this_param,min_inds);
end

end