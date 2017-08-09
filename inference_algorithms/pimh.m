function [samples, log_Zs, b_accept, mus] = pimh(sampling_functions,weighting_functions,...
                                                    N,resample_method,n_iter,b_compress,f,b_Rao_Black)
%pimh    Particle independent Metropolis Hastings
%
% Performs PIMH inference.  See section 4 in the iPMCMC paper or Particle
% Markove chain Monte Carlo methods, Andrieu et al (2010)
%
% Required inputs:
%   sampling_functions = See infer.m
%   weighting_functions = See infer.m
%   N (+ve integer) = Number of particles, also N in paper
%   resample_method = Method used in resampling.  See resample_particles.m.
%                     If empty takes default from resample_particles.m
%   n_iter (+ve integer) = Number of iterations
%   b_compress (boolean) = Whether to use compress_samples
%   f = Function to take expectation of.  Takes the var field of samples as
%       inputs.  See function_expectation.m.
%                               Default = []; (i.e. no estimate made)
%   b_Rao_Black (boolean) = Whether to Rao-Blackwellize and return all
%                           generated samples or just the retained particle
%
% Outputs:
%   samples = Object array of type stack_object containing details about
%             sampled variables, their weights and any constant variables
%   log_Zs = Marginal likelihood of individual sweeps
%   b_accept = Boolean vector indicating if that iteration is accepted
%   mus = Mean estimates of individual sweeps
%
% Tom Rainforth 08/06/16

log_Zs = NaN(n_iter,1);

if b_Rao_Black
    prop_sub_sample = 1;
else
    % Sampling the retained particle is equivalent to taking a 1/N proporiton 
    prop_sub_sample = 1/N;
end

b_compress = b_compress && b_Rao_Black;

[samples, log_Zs(1), mus] = smc_sweep(sampling_functions,weighting_functions,N,resample_method,b_compress,f,prop_sub_sample);
mus = repmat(mus,n_iter,1);

if ~b_compress
    % Memory management once have information from the first iteration
    [samples,b_compress] = memory_check(samples,n_iter,numel(sampling_functions));
end

samples = repmat(samples,n_iter,1);

b_accept = true(n_iter,1);

for iter=2:n_iter
    % Independent proposal
    [sample_proposed, log_Z_proposed, mu_proposed] = smc_sweep(sampling_functions,weighting_functions,...
                                                                N,resample_method,b_compress,f,prop_sub_sample);
    % MH acceptance ratio for new particle set
    bKeep = rand<min(1,exp(log_Z_proposed-log_Zs(iter-1)));
    if bKeep
        samples(iter) = sample_proposed;
        log_Zs(iter) = log_Z_proposed;
        mus(iter,:) = mu_proposed;
    else
        samples(iter) = samples(iter-1);
        log_Zs(iter) = log_Zs(iter-1);
        b_accept(iter) = false;
        mus(iter,:) = mus(iter-1,:);
    end
    
    if ~b_Rao_Black
        samples(iter).relative_particle_weights = 1;
    end
end

end