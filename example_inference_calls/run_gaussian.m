% Gaussian model with pgibbs inference, uncompressed output

% Do inference on the Gaussian model using particle gibbs
samples_g = infer('gaussian_model',[],'pgibbs','n_particles',1e4,'n_iter',1e3,'b_compress',false);

% Calculate mean and standard deviation
[mu_g_out,sig_g_out] = empirical_moments(samples_g,[1,2],'mu');

% Precalculated ground truths
mu_g_truth = 7.25;
sig_g_truth = sqrt(1/1.2);

gaussian_KL = log(sig_g_out)-log(sig_g_truth)+((sig_g_truth)^2+(mu_g_truth-mu_g_out)^2)/(2*(sig_g_out^2))-0.5;

disp(['Gaussian mu truth ' num2str(mu_g_truth) ' sampled ' num2str(mu_g_out)]);
disp(['Gaussian sig truth ' num2str(sig_g_truth) ' sampled ' num2str(sig_g_out)]);
disp(['Gaussian KL divergence to ground truth ' num2str(gaussian_KL)]);

% Plot output histogram
histogram_plotter(samples_g,'mu',100,1,1);

drawnow;