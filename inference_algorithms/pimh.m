function [samples, acceptance_ratio, log_Zs, b_accept] = pimh(sampling_functions,weighting_functions,...
                                                                N,n_iter,b_compress,b_sparse)
%pimh    Particle independent Metropolis Hastings
%
% Performs PIMH inference.  See section 4 in the iPMCMC paper or Particle
% Markove chain Monte Carlo methods, Andrieu et al (2010)
%
% Required inputs:
%   sampling_functions = See infer.m
%   weighting_functions = See infer.m
%   N (+ve integer) = Number of particles, also N in paper
%   n_iter = Number of independent sweeps to perform
%   b_compress (boolean) = Whether to use compress_samples
%
% Optional inputs:
%   b_sparse (boolean) = Takes the sparsity exploitation to the next level
%                        by compressing on-the-fly after each step in the
%                        state sequence.  Note that this incurs significant
%                        computational overhead and should be avoided
%                        whenever possible.
%                               Default = false;
%
% Outputs:
%   samples = Object array of type stack_object containing details about
%             sampled variables, their weights and any constant variables
%   acceptance_ratio = Average number of times the new sweep is accepted
%   log_Zs = Marginal likelihood of individual sweeps
%   b_accept = Boolean vector indicating if that iteration is accepted
%
% Tom Rainforth 08/06/16

if ~exist('b_sparse','var')
    b_sparse = false;
end

log_Zs = NaN(n_iter,1);

[samples, log_Zs(1)] = smc_sweep(sampling_functions,weighting_functions,N,b_compress,b_sparse);

if ~b_compress
    %% Memory management once have information from the first iteration
    S = whos('samples');
    s_mem = S.bytes*n_iter;
    if s_mem>5e7
        try
            memory_stats = memory;
            largest_array = memory_stats.MaxPossibleArrayBytes;
        catch
            % memory function is only availible in windows
            largest_array = 4e9;
        end
        
        if S.bytes*n_iter > (largest_array/20)
            warning('In danger of swamping memory and crashing, turning b_compress on');
            b_compress = true;
            samples = compress_samples(samples, numel(sampling_functions));
        end
    end
end

samples = repmat(samples,n_iter,1);

b_accept = true(n_iter,1);

for iter=2:n_iter    
    [sample_proposed, log_Z_proposed] = smc_sweep(sampling_functions,weighting_functions,N,b_compress,b_sparse);
    bKeep = rand<min(1,exp(log_Z_proposed-log_Zs(iter-1)));
    if bKeep
        samples(iter) = sample_proposed;
        log_Zs(iter) = log_Z_proposed;
    else
        samples(iter) = samples(iter-1);
        log_Zs(iter) = log_Zs(iter-1);
        b_accept(iter) = false;
    end
end

acceptance_ratio = mean(double(b_accept));

end