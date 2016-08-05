clear all
close all

ground_truths = load(['..' filesep() 'example_models' filesep() 'kalman_filter_data' filesep() 'ground_truth_summary.mat'],'truths');

for n=1:10;

data_lss = load(['..' filesep() 'example_models' filesep() 'kalman_filter_data' filesep() 'kalman_filter_data_' num2str(n) '.mat']);
ground_truth(n) = ground_truths.truths.(['b' num2str(n)]);

n_iter = 10000;
samples_lss_ipmcmc(n) = infer('kalman',data_lss.model_inputs,'ipmcmc','n_particles',100,'n_iter',n_iter,'M',32,'P',16,'b_parallel',true,'b_compress',true);
samples_lss_mPG(n) = infer('kalman',data_lss.model_inputs,'independent_nodes','n_particles',100,'n_iter',n_iter,'Ms',[32,0,0],'b_parallel',true,'b_compress',true);

h1 = figure;
h2 = figure;

plot_lss_error(samples_lss_ipmcmc(n),ground_truth(n),h1,h2,'r');
plot_lss_error(samples_lss_mPG(n),ground_truth(n),h1,h2,'b');

drawnow;

end