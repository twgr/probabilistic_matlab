function plot_lss_error(samples_lss,ground_truth,h1,h2,colour)

[mu_lss,sig_lss,ske_lss,kur_lss] = samples_lss.empirical_moments_convergence([1,2,3,4],'x');
n1 = size(mu_lss,1);
n2 = 3;
n3 = size(mu_lss,2)/n2;

mu_lss = reshape(mu_lss,n1,n2,n3); 
sig_lss = reshape(sig_lss,n1,n2,n3); 
ske_lss = reshape(ske_lss,n1,n2,n3); 
kur_lss = reshape(kur_lss,n1,n2,n3);

for n=1:3
    mu_error(:,n,:) = bsxfun(@minus,mu_lss(:,n,:),reshape(ground_truth.(['x_' num2str(n)]).Mean,[1,1,size(mu_lss,3)]));
    sig_error(:,n,:) = bsxfun(@minus,sig_lss(:,n,:),reshape(ground_truth.(['x_' num2str(n)]).Std_dev,[1,1,size(mu_lss,3)]));
    ske_error(:,n,:) = bsxfun(@minus,ske_lss(:,n,:),reshape(ground_truth.(['x_' num2str(n)]).Skewness,[1,1,size(mu_lss,3)]));
    kur_error(:,n,:) = bsxfun(@minus,kur_lss(:,n,:),reshape(ground_truth.(['x_' num2str(n)]).Excess_kurtosis,[1,1,size(mu_lss,3)]));
end

mu_conv = squeeze(mean(mean(mu_error.^2,3),2));
sig_conv = squeeze(mean(mean(sig_error.^2,3),2));
ske_conv = squeeze(mean(mean(ske_error.^2,3),2));
kur_conv = squeeze(mean(mean(kur_error.^2,3),2));

mu_pos = squeeze(mean(mu_error(end,:,:).^2,2));
sig_pos = squeeze(mean(sig_error(end,:,:).^2,2));
ske_pos = squeeze(mean(ske_error(end,:,:).^2,2));
kur_pos = squeeze(mean(kur_error(end,:,:).^2,2));

figure(h1);
subplot(2,2,1); loglog(mu_conv,colour); hold on;title('Mean');
xlabel('MCMC Iteration'); ylabel('Mean Squared Error');
subplot(2,2,2); loglog(sig_conv,colour); hold on;title('Std Dev');
xlabel('MCMC Iteration'); ylabel('Mean Squared Error');
subplot(2,2,3); loglog(ske_conv,colour); hold on;title('Skewness');
xlabel('MCMC Iteration'); ylabel('Mean Squared Error');
subplot(2,2,4); loglog(kur_conv,colour); hold on; title('Excess Kurtosis');
xlabel('MCMC Iteration'); ylabel('Mean Squared Error');

figure(h2);
subplot(2,2,1); semilogy(mu_pos,colour); hold on; title('Mean');
xlabel('State space time step t'); ylabel('Mean Squared Error');
subplot(2,2,2); semilogy(sig_pos,colour); hold on;  title('Std Dev');
xlabel('State space time step t'); ylabel('Mean Squared Error');
subplot(2,2,3); semilogy(ske_pos,colour); hold on; title('Skewness');
xlabel('State space time step t'); ylabel('Mean Squared Error');
subplot(2,2,4); semilogy(kur_pos,colour); hold on; title('Excess Kurtosis');
xlabel('State space time step t'); ylabel('Mean Squared Error');