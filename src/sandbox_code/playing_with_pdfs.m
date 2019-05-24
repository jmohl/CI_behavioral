%playing with probability distributions


%for unisensory condition
xa = -50:50;
sa(:,1) = -50:50;
A_sig = 5
prior_sig = 20;
S = [ -24 -18 -12 -6 6 12 18 24];
cond_a = [-24 -6 6 24];
prior_mu = 0;


%p(sa|xa)
likelihood = bsxfun_normpdf(xa,sa,A_sig);
%figure
%imagesc(likelihood);
% title('likelihood p(xa|sa)')

%pr(sa)
%prior = ones(size(Ta))/length(Ta);
prior = zeros(size(xa));
prior(ismember(xa,S)) = 1/length(S);
posterior = likelihood.*prior';
posterior = posterior./(1/sum(posterior,'all'));
%figure
%imagesc(posterior)
% xlabel('xa')
% ylabel('sa')
% title('p(sa|xa) posterior - discrete prior')

%p(sa) normal prior
prior_norm = normpdf(xa,0,prior_sig);
post_norm = bsxfun(@times,likelihood,prior_norm); %what I want this to do is multiply every row by the prior, but it's not really doing that.
post_norm = post_norm./sum(post_norm,2); %rescales every tow to make it a probability
%figure
%imagesc(post_norm);
% xlabel('xa')
% ylabel('sa')
% title('p(sa|xa) posterior - norm prior')

%mixture of normals prior
prior_norm2 = 0.5* normpdf(xa,15,prior_sig) + 0.5* normpdf(xa,-15,prior_sig);
post_norm2 = bsxfun(@times,likelihood,prior_norm2); %what I want this to do is multiply every row by the prior, but it's not really doing that.
post_norm2 = post_norm2./sum(post_norm2,2); %rescales every tow to make it a probability
%figure
%imagesc(post_norm2);
% xlabel('xa')
% ylabel('sa')
% title('p(sa|xa) posterior - double norm prior')



%in order to compare with real data, need to marginalize over xa values for
%a given Ta value
xpdf_A = bsxfun_normpdf(xa, cond_a', A_sig);
xpdf_A = bsxfun(@rdivide, xpdf_A, qtrapz(xpdf_A, 2)); %this only matters if the area of each slice is different than 1.
post_rescale(1,:,:) = posterior;
est_A = qtrapz(bsxfun(@times, xpdf_A, post_rescale), 2);
est_A = est_A./sum(est_A,3);
%figure
%imagesc(squeeze(est_A));
% title('p(ra|sa) discrete')

post_norm_rescale(1,:,:) = post_norm;
est_A_norm = qtrapz(bsxfun(@times, xpdf_A, post_norm_rescale), 2);
est_A_norm = est_A_norm./sum(est_A_norm,3);
%figure
%imagesc(squeeze(est_A_norm));
% title('p(ra|sa) normal')

post_norm2_rescale(1,:,:) = post_norm2;
est_A_norm2 = qtrapz(bsxfun(@times, xpdf_A, post_norm2_rescale), 2);
est_A_norm2 = est_A_norm2./sum(est_A_norm2,3);
%figure
%imagesc(squeeze(est_A_norm2));
% title('p(ra|sa) normal2')


figure
hold on
plot(xa,squeeze(est_A(1,:)))
plot(xa,squeeze(est_A(2,:)))
plot(xa,squeeze(est_A(3,:)))
plot(xa,squeeze(est_A(4,:)))
xticks(cond_a)
grid on
title('sensory response posteriors p(ra|sa), discrete prior')

figure
hold on
plot(xa,squeeze(est_A_norm(1,:)))
plot(xa,squeeze(est_A_norm(2,:)))
plot(xa,squeeze(est_A_norm(3,:)))
plot(xa,squeeze(est_A_norm(4,:)))
xticks(cond_a)
grid on
title('sensory response posteriors p(ra|sa), normal prior')

figure
hold on
plot(xa,squeeze(est_A_norm2(1,:)))
plot(xa,squeeze(est_A_norm2(2,:)))
plot(xa,squeeze(est_A_norm2(3,:)))
plot(xa,squeeze(est_A_norm2(4,:)))
xticks(cond_a)
grid on
title('sensory response posteriors p(ra|sa), dual normal prior')



%% one thing I really need to do is compare doing the prior + likelihood combination this way to solving it analytically.
As_mu = (xa/A_sig^2 + prior_mu/prior_sig^2) / (1/A_sig^2 + 1/prior_sig^2); %incorporating prior sig on unisensory estimates
As_sig = sqrt((1/A_sig^2 + 1/prior_sig^2)^-1);
A_pdf = bsxfun_normpdf(xa,As_mu',As_sig);
%figure
%imagesc(A_pdf)
A_pdf_rescale(1,:,:) = A_pdf;
est_A_pdf = qtrapz(bsxfun(@times, xpdf_A, A_pdf_rescale), 2);
%figure
%imagesc(squeeze(est_A_pdf))
%figure
%imagesc(squeeze(est_A_norm))

%figure
plot(xa,squeeze(est_A_pdf(2,:)/sum(est_A_pdf(2,:,:))))
hold on
plot(xa,squeeze(est_A_norm(2,:,:)/sum(est_A_norm(2,:,:))))


        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
