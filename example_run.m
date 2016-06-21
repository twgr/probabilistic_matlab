%% Branching with smc inference, compressed output

samples_branching = infer('branching_model',[],'smc','n_particles',5e5,'n_iter',3,'b_compress',true);

densities_branching_r = zeros(1,20);
densities_branching_r(1:max(samples_branching.var.r)+1) = results_binning(samples_branching,'r',1,size(samples_branching.var.r,2),0:1:max(samples_branching.var.r),true);
densities_branching_r_truth = [0.0209421460738810;0.120049724359027;0.0676937411893348;1.02017216903057e-09;1.00000000006526e-32;0.333292841871083;0.222305145210522;0.126917947983607;0.0633921064841567;0.0281583673328266;0.0112158251425169;0.00410248708336039;0.00135825297633444;0.000414215967343333;0.000116031661294933;3.05009697889260e-05;8.31844630607070e-06;1.49305446519218e-06;6.39880485082362e-07;2.13293495027454e-07];
figure;
bar(-0.125:1:18.875,densities_branching_r(1:20),0.25);
hold on
bar(0.125:1:19.125,densities_branching_r_truth(1:20),0.25,'r');
xlabel('r')
ylabel('P(r)')
title('Probabilities of r given by Branching model')
legend('Output','Truth');

%% Gaussian model with pgibbs inference, uncompressed output

samples_g = infer('gaussian_model',[],'pgibbs','n_particles',1e4,'n_iter',1e3,'b_compress',false);

mu_g_out = empirical_mean(samples_g,'mu');
mu_g_truth = 7.25;
sig_g_out = sqrt(empirical_covariance(samples_g,true,'mu'));
sig_g_truth = sqrt(1/1.2);
disp(['Gaussian mu truth ' num2str(mu_g_truth) ' sampled ' num2str(mu_g_out)]);
disp(['Gaussian sig truth ' num2str(sig_g_truth) ' sampled ' num2str(sig_g_out)]);

histogram_plotter(samples_g,'mu',1,1,1,100);

%% HMM with pimh inference

data_hmm = [0.9 0.8 0.7 0.0 -0.025 -5.0 -2.0 -0.1 0.0 0.13 0.45 6 0.2 0.3 -1 -1]';
samples_hmm = infer('hmm',data_hmm,'pimh','n_particles',1e4,'n_iter',20);

densities_hmm = results_binning(samples_hmm,'x',1,size(samples_hmm.var.x,2),[1 2 3],true);
densities_hmm_truth = [0.377510329270661,0.309531763819114,0.312957906910225;0.0414378861264862,0.406623078612501,0.551939035261013;0.0536324337803479,0.254680440661969,0.691687125557683;0.0456534970142020,0.230082297435319,0.724264205550479;0.105685915918161,0.120813663181357,0.773500420900481;0.0711054851992244,0.172144123412619,0.756750391388157;0.929464549092631,0.000123965467080903,0.0704114854402879;0.457535923680403,0.0456331420588863,0.496830934260711;0.0921655338249388,0.217034056521611,0.690800409653451;0.101544091004343,0.136235592179542,0.762220316816116;0.0986827336089095,0.157463653483870,0.743853612907220;0.178493790579279,0.219544259159816,0.601961950260906;6.25964028429687e-06,0.984915586803889,0.0150781535558268;0.112921375903980,0.167067325247519,0.720011298848501;0.0556364838062831,0.184851043327823,0.759512472865894;0.201502367616899,0.0473264197254615,0.751171212657639;0.254553838741815,0.0611172011550286,0.684328960103157];
densities_hmm_distance = sqrt(sum(sum((densities_hmm-densities_hmm_truth).^2))/(3*16));
disp(['HMM scaled distance to truth ' num2str(densities_hmm_distance)]);

%% Nonlinear state space model with dpg inference

data_nlss = load(['example_models' filesep() 'non_linear_state_space_data' filesep() 'Wn_equals_1' filesep() 'nonlinear_state_space_generate_data1.mat']);
n_steps = 20;
model_info.observations = data_nlss.Y(1:n_steps)';
n_iter = 20;
samples_nlss = infer('nonlinear_state_space',model_info,'ipmcmc','n_particles',10000,'n_iter',n_iter,'M',16,'P',5,'b_parallel',false);
 
histogram_plotter(samples_nlss,'x',5,4,1:20,300);
ESS = ess(samples_nlss,'x')/n_iter;
figure;
semilogy(ESS);
xlabel('Step in state space');
ylabel('Effective sample size per iteration');