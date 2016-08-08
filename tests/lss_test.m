clear all
close all

this_loc = fileparts(mfilename('fullpath'));

ground_truths = load([this_loc filesep() '..' filesep() 'example_models' filesep() 'kalman_filter_data' filesep() 'ground_truth_summary.mat'],'truths');

M = 32;
common_options = {'resample_method','multinomial','n_particles',100,'n_iter',1000,'b_parallel',true,'b_compress',true,'rng_seed',1};

for n=1:10;

data_lss = load([this_loc filesep() '..' filesep() 'example_models' filesep() 'kalman_filter_data' filesep() 'kalman_filter_data_' num2str(n) '.mat']);
ground_truth = ground_truths.truths.(['data_' num2str(n)]);

samples_lss_ipmcmc = infer('kalman',data_lss.model_inputs,'ipmcmc','M',M,common_options{:});
samples_lss_mPG = infer('kalman',data_lss.model_inputs,'independent_nodes','Ms',[M,0,0],common_options{:});

summary_lss_ipmcmc = samples_lss_ipmcmc.results_summary;
summary_lss_mPG = samples_lss_mPG.results_summary;

h1 = figure('units','normalized','outerposition',[0 0 1 1]);
h2 = figure('units','normalized','outerposition',[0 0 1 1]);

save(['lss_results_' num2str(n)], 'summary_lss_ipmcmc', 'summary_lss_mPG', 'ground_truth');

plot_lss_error(summary_lss_ipmcmc,ground_truth,h1,h2,'r');
plot_lss_error(summary_lss_mPG,ground_truth,h1,h2,'b');

drawnow;

end