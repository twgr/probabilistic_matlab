% Branching with smc inference, compressed output

% Do inference on the branching model using smc
samples_branching = ...
    infer('branching_model',[],'smc','n_particles',5e5,'n_iter',3,'b_parallel',true,'b_compress',true);

% Calculate probability for each value of r
p_r = zeros(1,20);
p_r(1:max(samples_branching.var.r)+1) = ...
    results_binning(samples_branching,'r',max(samples_branching.var.r)+1,true,[],[],0:1:max(samples_branching.var.r))';

% Precalculated ground truth
p_r_truth = [0.0209421460738810;0.120049724359027;0.0676937411893348;1.02017216903057e-09;1.00000000006526e-32;...
                               0.333292841871083;0.222305145210522;0.126917947983607;0.0633921064841567;0.0281583673328266;...
                               0.0112158251425169;0.00410248708336039;0.00135825297633444;0.000414215967343333;0.000116031661294933;...
                               3.05009697889260e-05;8.31844630607070e-06;1.49305446519218e-06;6.39880485082362e-07;2.13293495027454e-07];
                           
% Plot the results in a bar plot
figure('units','normalized','outerposition',[0 0 1 1]);
bar(-0.125:1:18.875,p_r(1:20),0.25);
hold on
bar(0.125:1:19.125,p_r_truth(1:20),0.25,'r');
xlabel('r')
ylabel('P(r)')
title('Probabilities of r given by Branching model')
legend('Output','Truth');
set(gca,'FontSize',32);

drawnow;