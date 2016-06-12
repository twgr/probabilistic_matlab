function [particles, log_Z, variables_step, sizes_step] = smc_sweep(sampling_functions,weighting_functions,N,b_compress,b_online_compression)
%smc_sweep   Carries out a single unconditional SMC sweep
%
% Performs sequential Monte Carlo (SMC) as per Algorithm 1 in the paper.
%
% Required inputs:
%   sampling_functions = See infer.m
%   weighting_functions = See infer.m
%   N (+ve integer) = Number of particles, also N in paper
%
% Optional inputs:
%   b_compress (boolean) = Whether to use compress_samples
%                               Default = false;
%   b_online_compression (boolean) = 
%                        Takes the sparsity exploitation to the next level
%                        by compressing on-the-fly after each step in the
%                        state sequence.  Note that this incurs significant
%                        computational overhead and should be avoided 
%                        whenever possible.
%                               Default = false;
%
% Outputs:
%   particles = Object of type stack_object storing all the samples.
%   log_Z = Log marginal likelihood estimate as per equation 4 in the paper
%   variables_step = Variables present at each step in state sequence
%   sizes_step = Sizes of variables at each step in state sequence
%
% Tom Rainforth 07/06/16

if ~exist('b_compress','var'); b_compress = false; end

if ~exist('b_online_compression','var'); b_online_compression = false; end

% Global parameter that allows for communication of the sample size with
% the sampling_functions and weighting_functions if desired.
global sample_size
sample_size = N;

particles = stack_object;
log_Z = 0;
[variables_step,sizes_step] = deal(cell(numel(sampling_functions),1));

if b_online_compression % Preallocation in sparse case    
    parents = sparse(N,numel(sampling_functions));
end

for n=1:numel(sampling_functions)
    % Simulate forwards as per eq 1a in the paper
    particles = sampling_functions{n}(particles);
    % Weight as per eq 1b in the paper
    log_weights = weighting_functions{n}(particles);
    
    % Store intermediary details
    variables_step{n} = fields(particles.var);
    for v=1:numel(variables_step{n})
        if isnumeric(particles.var.(variables_step{n}{v})) || (size(particles.var.(variables_step{n}{v}),2)>1)
            % Here we have an array so just need to store the second
            % dimension which is constant for all samples
            sizes_step{n}{v} = size(particles.var.(variables_step{n}{v}),2);
        elseif iscell(particles.var.(variables_step{n}{v}))
            % Here we have a n_samplesx1 size cell array so the size
            % might be different within each cell and needs storing
            % seperately
            sizes_step{n}{v} = cellfun(@(x) size(x,2), particles.var.(variables_step{n}{v}));
        else
            sizes_step{n}{v} = 1;
        end
    end
    
    if n~=numel(sampling_functions)
        % When not at the final step, perform resampling.
        [i_resample, log_Z_step, n_times_sampled] = resample_step(log_weights,numel(log_weights));
        
        assert(~isnan(log_Z_step),'log_Z_step is NaN');
        log_Z = log_Z+log_Z_step;
              
        if ~b_online_compression
            % Update particles object as per the resampling indices
            p_fields = fields(particles.var);
            for n_f = 1:numel(p_fields)
                particles.var.(p_fields{n_f}) = particles.var.(p_fields{n_f})(i_resample,:);
            end
        else
%% Online compression code
% Compresses after each sample using a sparse array and implicitly storing
% the ancestral path in the structure of the sparse array.
            n_times_sampled = n_times_sampled(1:end-1);
            b_sampled = n_times_sampled~=0;
            n_times_sampled_no_zeros = n_times_sampled(b_sampled);
            heritage = find(b_sampled);
            sparse_fields = fields(particles.var);
            for n_f = 1:numel(sparse_fields)
                assert(isnumeric(particles.var.(sparse_fields{n_f})),'Sparse format currently only supports numeric variables');
                assert(size(particles.var.(sparse_fields{n_f}),2)==1,'In sparse format currently do not support variables being more than one dimensional');
                particles.sparse_history.(sparse_fields{n_f})(:,n) = sparse(heritage,1,particles.var.(sparse_fields{n_f})(heritage),N,1);
                particles.var.(sparse_fields{n_f}) = particles.var.(sparse_fields{n_f})(i_resample,:);
            end
            parents(:,n) = i_resample; %#ok<SPRIX>
            for n_a = n:-1:2
                if n_a~=1
                    n_times_sampled_local = accumarray(heritage, n_times_sampled_no_zeros, [N,1], [], [], 1);
                    i_sampled = find(n_times_sampled_local);
                    n_times_sampled_no_zeros = nonzeros(n_times_sampled_local);
                    heritage = full(parents(i_sampled,n_a-1));
                end
                
                i_sampled_reduced = parents(i_sampled,n_a-1);
                parents(:,n_a-1) = sparse(i_sampled,1,i_sampled_reduced,N,1); %#ok<SPRIX>
                for n_f = 1:numel(sparse_fields)
                    i_sampled_reduced = sort(i_sampled_reduced);
                    b_diff = [true;diff(i_sampled_reduced)~=0];
                    i_sampled_reduced = i_sampled_reduced(b_diff);
                    particles.sparse_history.(sparse_fields{n_f})(:,n_a-1) = ...
                        sparse(i_sampled_reduced,1,particles.sparse_history.(sparse_fields{n_f})(i_sampled_reduced,n_a-1),N,1);
                end
            end
        end
        
    end
end

% Set outputs
z_max = max(log_weights);
w = exp(log_weights-z_max);
log_Z = log_Z+z_max+log(sum(w))-log(numel(w)); % Marginal likelihood
particles.relative_particle_weights = w/sum(w);

%% Tying up of sample compression as required.

if b_online_compression
    parents(:,end) = sparse((1:N)');
    heritage = full(parents(:,end));
    weights_no_zeros = particles.relative_particle_weights;
    weights_no_zeros(weights_no_zeros==0) = min(1e-32,min(weights_no_zeros(weights_no_zeros~=0))/1e10);
    sparse_relative_weights = sparse(N,numel(sampling_functions));
    for n=numel(sampling_functions):-1:1
        weights_local = accumarray(heritage, weights_no_zeros, [N,1], [], [], 1);
        i_sampled = find(weights_local);
        weights_no_zeros = nonzeros(weights_local);
        sparse_relative_weights(:,n) = accumarray(i_sampled, weights_no_zeros, [N,1], [], [], 1); %#ok<SPRIX>
        if n~=1
            heritage = full(parents(i_sampled,n-1));
        end
    end
    
    for n_f = 1:numel(sparse_fields)
        particles.sparse_history.(sparse_fields{n_f})(:,numel(sampling_functions)) = sparse(particles.var.(sparse_fields{n_f}));
        particles.var.(sparse_fields{n_f}) = particles.sparse_history.(sparse_fields{n_f});
        particles.sparse_history = rmfield(particles.sparse_history,sparse_fields{n_f});
    end
end

if b_compress
    particles = compress_samples(particles,numel(weighting_functions));
end

if exist('sparse_relative_weights','var')
    particles.sparse_variable_relative_weights = sparse_relative_weights;
end

end
