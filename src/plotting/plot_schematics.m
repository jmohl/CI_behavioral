%% Schematic plots for CI and location estimation
%
% -------------------
% Jeff Mohl
% 7/9/19
% -------------------
%
% Description: Making schematic illustrations for comparison with data.
% plot1 demonstrate CI using probability distributions
% plot2 should be illustration of unity judgement under different
% possibilities
% plot3 should be illustration of localization plot under CI and non-CI
% models

%% plot 1: demonstration distributions for localization of vis and aud tars
figure
%options
A_sig = 7;
V_sig = 4;
A_tar = 10;
V_tar = -10;
bin_size = .1;
xrange = -50:bin_size:50;
%make pdfs
A_pdf = normpdf(xrange,A_tar,A_sig);
V_pdf = normpdf(xrange,V_tar,V_sig);
AV_pdf = A_pdf .* V_pdf; %optimal integration assuming normal dists and no prior
AV_pdf = AV_pdf ./ (sum(AV_pdf)*bin_size); %normalize to sum to 1.
%make plot
plot(xrange,A_pdf,'linewidth',2,'color','r');
hold on
plot(xrange,V_pdf,'linewidth',2,'color','b');
plot(xrange,AV_pdf,'linewidth',2,'color','k');
legend('Auditory', 'Visual', 'Integrated');
xlabel('Sensory parameter (location)');
ylabel('Probability');

% now make that plot, but show different lines for what would be expected
% under causal inference, i.e. different lines for segregated, integrated 

%% plot 1 again, but this time using the code I've written for the actual model

%options
A_sig = 5;
V_sig = 3;
A_tar = [6;12];
V_tar = [-6;-12];
bin_size = .5;
prior_mu = 0;
prior_sig =20;
p_common = .8;
xrange = -50:bin_size:50;
xrange_A(1,1,:) = xrange;
xrange_V(1,:,1) = xrange;

%all of this is just taken from the datalike code, might be simpler to just
%run it once and get the output that way, but it expects data. 
c1post = get_c1post(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,p_common); %switched prior mu with 0 here (3rd param) because I'm being lazy, will remove once I add the new CI priors.
int_pdf = get_integrate_pdf(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,xrange);
[~,A_seg_pdf,V_seg_pdf] = get_segregate_pdf(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,xrange); %seg pdf is 1x1x(xA)x(eval_range) array, so pdf is in 4th dimension for a given value of xA

AV_c1_pdf = bsxfun(@times,c1post, int_pdf) ; %splitting these because my responses are split between 1 and 2 saccade cases.
V_c2_pdf = bsxfun(@times,(1-c1post),V_seg_pdf);
A_c2_pdf = bsxfun(@times,(1-c1post), A_seg_pdf);

xpdf_V = bsxfun_normpdf(xrange_V, V_tar,V_sig);
xpdf_V = bsxfun(@rdivide, xpdf_V, qtrapz(xpdf_V, 2)); % Not multiplying by volume element, this doesn't matter if the grid volumn = 1
xpdf_A = bsxfun_normpdf(xrange_A, A_tar, A_sig);
xpdf_A = bsxfun(@rdivide, xpdf_A, qtrapz(xpdf_A, 3));  % Not multiplying by volume element

prmat_AV_c1 = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), AV_c1_pdf), 2), 3); %this is a recreation of what VestBMS_finalqtrapz does but not done in C.
prmat_A_c2 = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), A_c2_pdf), 2), 3);

prmat_V_c2 = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), V_c2_pdf), 2), 3);
prmat_V_shape(:,:,1) = prmat_V_c2;
prmat_A_shape(:,1,:) = prmat_A_c2;
prmat_AV_c2 = bsxfun(@times,prmat_A_shape,prmat_V_shape);%problem here is normalizing to make the probabilities match, because this is going to significantly reduce them
%renormalize so that total probability will sum to 1 for each condition
%after adding c=1 case
prmat_AV_c2 = bsxfun(@rdivide,prmat_AV_c2,sum(prmat_V_c2,4));

%diagonalize
diag_array = diag(ones(length(xrange),1));
diag_array = logical(diag_array);
prmat_sac = prmat_AV_c2;
%add C=1 case to diagonal, producing final cxSvxSa matrix
prmat_sac(:,diag_array) = prmat_sac(:,diag_array) + squeeze(prmat_AV_c1(:,1,1,:));

%sum over A or V locs, to get A or V pdf in 2d
figure
 set(gcf,'Position',[25,100,1600,500])
for plot_pair = [1,2]
    subplot(1,2,plot_pair);
    hold on

cA = squeeze(sum(prmat_sac(plot_pair,:,:),2));
cV = squeeze(sum(prmat_sac(plot_pair,:,:),3));
cAV = squeeze(prmat_sac(plot_pair,diag_array));
%if want to plot single saccades separately, then include this
% cA = cA' - cAV;
% cV = cV - cAV;
%normalize all pdfs to 1, this eliminates the ratio but makes it easier to
%compare variance and mean of all
cA = cA/sum(cA);
cV=cV/sum(cV);
cAV = cAV/sum(cAV);

plot(xrange,cA,'r')
plot(xrange,cV,'b')
% plot(xrange,cAV,'k')

[mv_V,mi_V]= max(cV);
[mv_A,mi_A]= max(cA);

% get fully segregated estimates
prmat_A_seg = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), A_seg_pdf), 2), 3);
prmat_V_seg = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), V_seg_pdf), 2), 3);

cA_seg= squeeze(prmat_A_seg(plot_pair,:,:,:));
cV_seg= squeeze(prmat_V_seg(plot_pair,:,:,:));
cA_seg=cA_seg/sum(cA_seg);
cV_seg=cV_seg/sum(cV_seg);

plot(xrange,cA_seg,'r--')
plot(xrange,cV_seg,'b--')
title(sprintf('%d aud, %d vis',A_tar(plot_pair),V_tar(plot_pair)));
ylabel('Probability')
xlabel('Stimulus feature (location)')
legend('Auditory CI','Visual CI','Integrated','Auditory unimodal','Visual Unimodal')
end





