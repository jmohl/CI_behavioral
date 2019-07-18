%% Schematic plots for CI and location estimation
%
% -------------------
% Jeff Mohl
% 7/9/19
% -------------------
%
% Description: Making schematic illustrations for comparison with data.
% I decided to make these plots by making the same PDFs I use in the model
% fitting process and making the plots out of those. This has the advantage
% of being a true reflection of what my code does, but there is kindof a
% disadvantage in that everything is binned and so the schematics are not
% going to have exact values. But I don't think that matters. It's also a
% little slow because it runs several integrations but I don't think that
% matters either.

%plot 1: pdfs for example conditions. Currently pdfs are normalized for
%comparison purposes and the single saccade trials (diagonal) are
%subtracted from the double saccade trials. You can remove this subtraction
%piece and end up with the kind of multi-modal distribution that is
%representative of paradigms which dont share this dual report scheme.

%plot 2: show unimodal and integrated bias, collapsing across conditions
%and plotting by target separation. This is the best way I've come up with
%so far to show the effects on localization in a single plot.

%plot 3: unity judgement by target sep. This plot probably won't make the
%paper because it conveys very little info but might be useful for
%demonstration purposes.
% clear;close all;
% 
% local_directory = 'C:\Users\jtm47\Documents\Projects\CI_behavioral\';
% cd(local_directory)
% addpath('src', 'src\plotting','results','data');
% figpath = 'results\figures\schematics';
% try
%     mkdir(figpath)
% end
function plot_schematics(figpath)
%% options
combination_rule = 1; %1 = bayes reweight, 2 = model selection 
A_sig = 5;
V_sig = 3;
A_tar = [-15;-12;-9;-6;-3;0;3;6;9;12;15];
V_tar = [15;12;9;6;3;0;-3;-6;-9;-12;-15];
bin_size = .25;
prior_mu = 0;
prior_sig =200; %easier to just make this very large rather than take out of code
p_common = .5;
xrange = -50:bin_size:50;
%% Generate pdfs
xrange_A(1,1,:) = xrange;
xrange_V(1,:,1) = xrange;

%all of this is just taken from the datalike code, might be simpler to just
%run it once and get the output that way, but it expects data.
c1post = get_c1post(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,p_common); %switched prior mu with 0 here (3rd param) because I'm being lazy, will remove once I add the new CI priors.
int_pdf = get_integrate_pdf(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,xrange);
[~,A_seg_pdf,V_seg_pdf] = get_segregate_pdf(xrange_A,xrange_V,prior_mu,A_sig,V_sig,prior_sig,xrange); %seg pdf is 1x1x(xA)x(eval_range) array, so pdf is in 4th dimension for a given value of xA

%choose weighting of pdfs based on CI strategy, reweighting or model
%selection
switch combination_rule
    case 1
        w_unity = c1post;
    case 2
        w_unity = zeros(size(c1post));
        w_unity(c1post > 0.5) = 1;
end

AV_c1_pdf = bsxfun(@times,w_unity, int_pdf) ;
V_c2_pdf = bsxfun(@times,(1-w_unity),V_seg_pdf);
A_c2_pdf = bsxfun(@times,(1-w_unity), A_seg_pdf);

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

% get fully segregated estimates
prmat_A_seg = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), A_seg_pdf), 2), 3);
prmat_V_seg = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), V_seg_pdf), 2), 3);

%% schematic 1, showing pdfs for 2 conditions (including bias)
figure
set(gcf,'Position',[25,100,1600,500])
plot_ind = 1;
for plot_pair = [8,10]
    subplot(1,2,plot_ind);
    hold on
    cA = squeeze(sum(prmat_sac(plot_pair,:,:),2));
    cV = squeeze(sum(prmat_sac(plot_pair,:,:),3));
    cAV = squeeze(prmat_sac(plot_pair,diag_array));
    %if want to plot single saccades separately, then include this
    cA = cA' - cAV;
    cV = cV - cAV;
    %normalize all pdfs to 1, this eliminates the ratio but makes it easier to
    %compare variance and mean of all
%     cA = cA/sum(cA);
%     cV=cV/sum(cV);
%     cAV = cAV/sum(cAV);
%     
    plot(xrange,cAV,'k')
    plot(xrange,cA,'r')
    plot(xrange,cV,'b')
    
    cA_seg= squeeze(prmat_A_seg(plot_pair,:,:,:));
    cV_seg= squeeze(prmat_V_seg(plot_pair,:,:,:));
%     cA_seg=cA_seg/sum(cA_seg);
%     cV_seg=cV_seg/sum(cV_seg);
    
    plot(xrange,cA_seg,'--','Color',[1,0,0,.25])
    plot(xrange,cV_seg,'--','Color',[0,0,1,.25])
    title(sprintf('%d aud, %d vis',A_tar(plot_pair),V_tar(plot_pair)));
    ylabel('Probability')
    xlabel('Stimulus feature (location)')
    legend('Integrated estimate (1 cause)','Auditory estimate (2 cause)','Visual estimate (2 cause)','Auditory Unimodal','Visual Unimodal')
    plot_ind = plot_ind+1;
end
saveas(gcf,sprintf('%s\\ex_pdfs',figpath),'svg');
%% schematic 2, format for collapsing across target conditions, showing only the mean of the distributions.
% there are two ways to do this. One is I can estimate the mean using the
% pdfs, but this does not really actually get the mean. The other is to
% calculate the means of the distributions, but this doesn't really work
% for the synthetic plots I'm making right now, because they dont have any
% data obviously.

% this also brings up a problem I'm going to have with the actual plots.
% Plotting the mean bias for various conditions doesn't really get at the
% heart of the issue, because the distributions are not normal. That's like
% kindof the point. Except that once you split the single and double
% saccade distributions this is much less of a problem, so it's probably
% fine.
cA = squeeze(sum(prmat_sac(:,:,:),2));
cV = squeeze(sum(prmat_sac(:,:,:),3));
cAV = squeeze(prmat_sac(:,diag_array));
cA_seg= squeeze(prmat_A_seg(:,:,:,:));
cV_seg= squeeze(prmat_V_seg(:,:,:,:));
%separate single and dual saccades
cA = cA - cAV;
cV = cV - cAV;

%interpolation code, slightly more accurate than using the max
respV = sum(cV.*xrange,2)./sum(cV,2);
respA = sum(cA.*xrange,2)./sum(cA,2);
respAV = sum(cAV.*xrange,2)./sum(cAV,2);
respA_seg = sum(cA_seg.*xrange,2)./sum(cA_seg,2);


AV_dif = A_tar - V_tar;
A_bias = respA - A_tar;
V_bias = respV - V_tar ;
AV_bias_A = respAV - A_tar;
AV_bias_V = respAV - V_tar ;
A_bias_seg = respA_seg -  A_tar;

figure
plot(AV_dif,A_bias,'.r-')
hold on
plot(AV_dif,V_bias,'.b-')
plot(AV_dif,AV_bias_A,'.k-')
% plot(AV_dif,AV_bias_V,'.k-')
plot(AV_dif,A_bias_seg,'.g--')
xlabel({'AV (A - V)';'<<<< A target to left of V | A target to right of V >>>>'})
ylabel({'Bias from target location';'(resp - target)'})
legend('A Saccade','V saccade','Single Saccade relative to A tar','A and V, no interaction','location', 'best')
title('Response bias vs target separation')
saveas(gcf,sprintf('%s\\means_by_sep',figpath),'svg');


%% schematic 3, unity judgement case
%this one is the most intuitive and therefore probably doesn't need and
%explanation as much. Still It will be good to have for a conversation with
%jenni so I can kindof lay out my ideas for the paper.

%the other thing that is a little annoying about this figure is that there
%really isn't much to it, and the difference between the models are super
%subtle with respect to this one (except the null model) so I can't even
%really use it to explain much. So probably won't include but whatever.
% 
% figure
% prmat_unity = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_V,xpdf_A), c1post), 2), 3); %this is a recreation of what VestBMS_finalqtrapz does but not done in C.
% plot(AV_dif, prmat_unity,'k.-')
% title('Percent report unity vs target separation')
% xlabel('\Delta AV (A - V)')
% ylabel('Percent report unity')
% saveas(gcf,sprintf('%s\\unity',figpath),'svg');
% 
% 

end

