% Nonlinear state space model with ipmcmc inference from the paper

dataset_to_use = 1;  % Scalar between 1 and 10

% Load data
data_nlss = load(['example_models' filesep() 'non_linear_state_space_data' filesep() 'nonlinear_state_space_generate_data' num2str(dataset_to_use) '.mat']);

% Choose number of steps
n_steps = 20;
model_info.observations = data_nlss.Y(1:n_steps)';

% Choose number of iterations
n_iter = 20;

% Run inference
samples_nlss = infer('nonlinear_state_space',model_info,'ipmcmc','n_particles',10000,'n_iter',n_iter,'M',16,'P',8,'b_parallel',true,'b_compress',true);
 
% Plot output histograms
histogram_plotter(samples_nlss,'x',300,4,5,false,1:20);

drawnow;

% Plot change in effective sample size over position in the state sequence
ESS = ess(samples_nlss,'x')/n_iter;
figure('units','normalized','outerposition',[0 0 1 1]);
semilogy(ESS);
xlabel('Step in state space');
ylabel('Effective sample size per iteration');
title('Effective sample size for NLSS Model');
set(gca,'FontSize',32);

drawnow;