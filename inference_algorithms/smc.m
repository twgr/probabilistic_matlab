function [samples, log_Z_total, log_Zs] = smc(sampling_functions,weighting_functions,N,n_iter,b_compress,b_sparse)
%smc    Independent SMC algorithm
%
% Performs SMC inference using independent sweeps.  See section 2.1 in the
% paper
%
% Required inputs:
%   sampling_functions = See infer.m
%   weighting_functions = See infer.m
%   N (+ve integer) = Number of particles, also N in paper
%   n_iter = Number of independent sweeps to perform
%
% Optional inputs:
%   b_compress (boolean) = Exploits the degeneracy caused by resampling to
%                          store the output using sparse matrices and an
%                          expliticly stored, sparse, ancestral lineage.
%                               Default = false;
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
%   log_Z_total = Marginal likehood of all sweeps combined
%   log_Zs = Marginal likelihood of individual sweeps
%
% Tom Rainforth 07/06/16

if b_sparse && b_compress
    warning('b_sparse does its own compression, no extra compression required');
    b_compress = false;
end

log_Zs = NaN(n_iter,1);

for iter=1:n_iter
    [samples(iter), log_Zs(iter)] = smc_sweep(sampling_functions,weighting_functions,N,b_compress,b_sparse); %#ok<AGROW>
    
    if iter==1
        %% Memory management once have information from the first iteration
        if ~b_compress

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
    end
end

%% Sort out the weights as per the marginal likelihoods
max_log_Z = max(log_Zs);
w = exp(log_Zs-max_log_Z);
log_Z_total = log(mean(w))+max_log_Z;
w = w/sum(w);

if ~b_sparse
    for iter=1:n_iter
        samples(iter).relative_particle_weights = w(iter)*samples(iter).relative_particle_weights;
    end
end

%% When sparse need to update the sparse weights
if ~isempty(samples(1).sparse_variable_relative_weights)
    if isstruct(samples(1).sparse_variable_relative_weights)
        sparse_fields = fields(samples(1).sparse_variable_relative_weights);
        for iter=1:n_iter
            for n_f=1:numel(sparse_fields)
                samples(iter).sparse_variable_relative_weights.(sparse_fields{n_f}) = samples(iter).sparse_variable_relative_weights.(sparse_fields{n_f})*w(iter);
            end
        end
    else
        for iter=1:n_iter
            samples(iter).sparse_variable_relative_weights = samples(iter).sparse_variable_relative_weights*w(iter);
        end
    end
end

end