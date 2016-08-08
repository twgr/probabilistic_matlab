% Linear Gaussian state space model from the paper with 32 independent
% particle gibbs inferences and iPMCMC with M=32 and P=16;

dataset_to_use = 1;

this_loc = fileparts(mfilename('fullpath'));

ground_truths = load([this_loc filesep() '..' filesep() 'example_models' filesep() 'kalman_filter_data' filesep() 'ground_truth_summary.mat'],'truths');
data_lss = load([this_loc filesep() '..' filesep() 'example_models' filesep() 'kalman_filter_data' filesep() 'kalman_filter_data_' num2str(dataset_to_use) '.mat']);
ground_truth = ground_truths.truths.(['data_' num2str(dataset_to_use)]);

samples_lss_ipmcmc = infer('kalman',data_lss.model_inputs,'ipmcmc','M',32,'P',16,'n_particles',100,'n_iter',200,'b_parallel',true,'b_compress',true);
samples_lss_mPG = infer('kalman',data_lss.model_inputs,'independent_nodes','Ms',[32,0,0],'n_particles',100,'n_iter',200,'b_parallel',true,'b_compress',true);

summary_lss_ipmcmc = samples_lss_ipmcmc.results_summary;
summary_lss_mPG = samples_lss_mPG.results_summary;

h1 = figure('units','normalized','outerposition',[0 0 1 1]);
h2 = figure('units','normalized','outerposition',[0 0 1 1]);

plot_lss_error(summary_lss_ipmcmc,ground_truth,h1,h2,'r');
plot_lss_error(summary_lss_mPG,ground_truth,h1,h2,'b');

drawnow;