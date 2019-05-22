%playing with probability distributions


%for unisensory condition
xa = -50:50;
A_sig = 5;
prior_sig = 3;
Ta = [-24 -6 6 24];
prior_mu = 0;


%p(sa|xa)
likelihood = bsxfun_normpdf(xa,xa',A_sig);
figure
imagesc(likelihood);

%pr(sa)
%prior = ones(size(Ta))/length(Ta);
prior = zeros(size(xa));
prior(ismember(xa,Ta)) = 1/length(Ta);
figure
posterior = likelihood.*prior';
imagesc(posterior)

%p(sa) normal prior
prior_norm = normpdf(xa,0,prior_sig);
post_norm = likelihood.*prior_norm;
figure
imagesc(post_norm);

%in order to compare with real data, need to marginalize over xa values for
%a given Ta value
xpdf_A = bsxfun_normpdf(xa, Ta', A_sig);
xpdf_A = bsxfun(@rdivide, xpdf_A, qtrapz(xpdf_A, 2)); %this only matters if the area of each slice is different than 1.
post_rescale(1,:,:) = posterior;
est_A = qtrapz(bsxfun(@times, xpdf_A, post_rescale), 2);
figure
imagesc(squeeze(est_A));

post_norm_rescale(1,:,:) = post_norm;
est_A_norm = qtrapz(bsxfun(@times, xpdf_A, post_norm_rescale), 2);
figure
imagesc(squeeze(est_A_norm));

figure
plot(xa,squeeze(est_A(1,:)))
hold on
plot(xa,squeeze(est_A_norm(1,:)))



%% one thing I really need to do is compare doing the prior + likelihood combination this way to solving it analytically.
As_mu = (xa/A_sig^2 + prior_mu/prior_sig^2) / (1/A_sig^2 + 1/prior_sig^2); %incorporating prior sig on unisensory estimates
As_sig = sqrt((1/A_sig^2 + 1/prior_sig^2)^-1);
A_pdf = bsxfun_normpdf(xa,As_mu',As_sig);
figure
imagesc(A_pdf)
A_pdf_rescale(1,:,:) = A_pdf;
est_A_pdf = qtrapz(bsxfun(@times, xpdf_A, A_pdf_rescale), 2);
figure
imagesc(squeeze(est_A_pdf))
figure
imagesc(squeeze(est_A_norm))

figure
plot(xa,squeeze(est_A_pdf(2,:)/sum(est_A_pdf(2,:,:))))
hold on
plot(xa,squeeze(est_A_norm(2,:,:)/sum(est_A_norm(2,:,:))))


        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
