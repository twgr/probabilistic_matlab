function [samples, log_Zs] = smc(sampling_functions,weighting_functions,...
                        N,n_iter,b_compress,b_parallel,n_islands,b_online_compression)
%smc    Independent SMC algorithm
%
% Performs SMC inference using independent sweeps.  See section 2.1 in the
% paper.  Independent smc runs produce consistent and unbiased samples.
% However, using evenly weighted "islands" (supported by setting
% n_islands>1) gives a lower variance but biased estimate.  This can
% therefore not be used in, for example, a nested SMC scheme.  Different 
% columns in the output will represent the different nodes (though if using
% infer these are later collapsed by process_final_samples if using
% infer.m).  Estimates are consistent as N->infinity or n_iter->infinity
% but not as n_islands->infinity.
%
% Inputs:
%   sampling_functions = See infer.m
%   weighting_functions = See infer.m
%   N (+ve integer) = Number of particles, also N in paper
%   n_iter = Number of independent sweeps to perform
%   b_compress (boolean) = Whether to use compress_samples
%   b_parallel (boolean) = Whether to run sweeps in parallel (number of
%                          threads automatically set by parfor)
%   n_islands (+ve integer) = Number of equally weighted islands to use as
%                        explained above.  The standard independent smc
%                        algorithm corresponds to n_islands=1 (Default = 1)
%   b_online_compression (boolean) = 
%                        Takes the sparsity exploitation to the next level
%                        by compressing on-the-fly after each step in the
%                        state sequence.  Note that this incurs significant
%                        computational overhead and should be avoided
%                        whenever possible. (Default = false)
%
% Outputs:
%   samples = Object array of type stack_object containing details about
%             sampled variables, their weights and any constant variables
%   log_Zs = Marginal likelihood of individual sweeps
%
% Tom Rainforth 12/06/16

if ~exist('n_islands','var')
    n_islands = 1;
end

if ~exist('b_online_compression','var')
    b_online_compression = false;
end

if b_online_compression && b_compress
    warning('b_online_compression does its own compression, no extra compression required');
    b_compress = false;
end

n_total = n_iter*n_islands;

% First sweep to use in initialization
log_Zs = NaN(n_total,1);
[samples, log_Zs(1)] = smc_sweep(sampling_functions,weighting_functions,N,b_compress,b_online_compression);

% Memory management
if ~b_compress    
    S = whos('samples');
    s_mem = S.bytes*n_total;
    if s_mem>5e7
        try
            memory_stats = memory;
            largest_array = memory_stats.MaxPossibleArrayBytes;
        catch
            % memory function is only availible in windows
            largest_array = 4e9;
        end
        
        if S.bytes*n_total > (largest_array/20)
            warning('In danger of swamping memory and crashing, turning b_compress on');
            b_compress = true;
            samples = compress_samples(samples, numel(sampling_functions));
        end
    end
end
samples = repmat(samples,n_total,1);

% Carry out independent sweeps
if b_parallel
    parfor iter=2:n_total
        [samples(iter), log_Zs(iter)] = smc_sweep(sampling_functions,weighting_functions,N,b_compress,b_online_compression);
    end
else
    for iter=2:n_total
        [samples(iter), log_Zs(iter)] = smc_sweep(sampling_functions,weighting_functions,N,b_compress,b_online_compression);
    end
end

%% Calculate relative weights and adjust field respectively

% Split into seperate islands
samples = reshape(samples,n_iter,n_islands);
log_Zs = reshape(log_Zs,n_iter,n_islands);
max_log_Z = max(log_Zs,[],1);
w = exp(bsxfun(@minus,log_Zs,max_log_Z));
% Relative weights of each island taken to be equal
w = bsxfun(@rdivide,w,sum(w,1))*(1/n_islands);

if ~isempty(samples(1).relative_particle_weights)
    for iter=1:n_total
        samples(iter).relative_particle_weights = w(iter)*samples(iter).relative_particle_weights;
    end
end

% When sparse also need to update the sparse weights
if ~isempty(samples(1).sparse_variable_relative_weights)
    if isstruct(samples(1).sparse_variable_relative_weights)
        sparse_fields = fields(samples(1).sparse_variable_relative_weights);
        for iter=1:n_total
            for n_f=1:numel(sparse_fields)
                samples(iter).sparse_variable_relative_weights.(sparse_fields{n_f}) = samples(iter).sparse_variable_relative_weights.(sparse_fields{n_f})*w(iter);
            end
        end
    else
        for iter=1:n_total
            samples(iter).sparse_variable_relative_weights = samples(iter).sparse_variable_relative_weights*w(iter);
        end
    end
end

end