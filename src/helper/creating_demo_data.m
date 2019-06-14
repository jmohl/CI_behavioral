%% create demo pdf for making schematic plots
%
% -------------------
% Jeff Mohl
% 6/14/19
% -------------------
%
% Description: generates model structures with likelihood values for a
% demonstration figures.

model_list = {[1 1 2 1]};
%setting fitting procedure options
MAXRNG = 50;
binsize = .25; %size of bins used for responses, in degrees.
fitoptions.make_plots = 0;
fitoptions.parameter_names = {'A_sig','V_sig','prior_sig','p_common','lambda_uni','lambda_loc','prior_mu'};
fitoptions.dynamic_bins = 0; %experimental

%parameter values
conditions = [-12,12;-12, 0];
responses = 1;
theta = [5,4,50,.5,.05,0,15];
eval_range = linspace(-MAXRNG,MAXRNG,MAXRNG*2/binsize + 1);
eval_midpoints = linspace(-MAXRNG+binsize/2,MAXRNG-binsize/2,length(eval_range)-1);

model = model_list{1};

[nll,prmat] = datalike(conditions,responses,theta,model,eval_midpoints);

this_dist = squeeze(prmat(1,:,:));
locations = eval_midpoints;
I_mat = logical(eye([length(locations),length(locations)]));
aud_dist = sum(this_dist,1)*binsize;
vis_dist = sum(this_dist,2)*binsize;
AV_dist = this_dist(I_mat)/sum(this_dist(I_mat),'all');
figure
hold on
plot(locations,aud_dist,'r','LineWidth',2)
plot(locations,vis_dist,'b','LineWidth',2)

plot(locations,AV_dist,'k','LineWidth',2)

figure