function [particles, log_Z, relative_weights] = smc_sweep(sampling_functions,weighting_functions,n_particles,b_compress,b_sparse)

if ~exist('b_sparse','var')
    b_sparse = false;
end

global sample_size
sample_size = n_particles;

particles = stack_object;
log_Z = 0;

if b_sparse    
    parents = sparse(n_particles,numel(sampling_functions));
end

for n=1:numel(sampling_functions)
    particles = sampling_functions{n}(particles);
    log_weights = weighting_functions{n}(particles);
    
    if n~=numel(sampling_functions)
        [i_resample, log_Z_step, n_times_sampled] = resample_step(log_weights,numel(log_weights));
           
        assert(~isnan(log_Z_step),'log_Z_step is NaN');
        log_Z = log_Z+log_Z_step;
        
        if b_sparse
            n_times_sampled = n_times_sampled(1:end-1);
            b_sampled = n_times_sampled~=0;
            n_times_sampled_no_zeros = n_times_sampled(b_sampled);
            heritage = find(b_sampled);
            sparse_fields = fields(particles.var);
            for n_f = 1:numel(sparse_fields)
                assert(isnumeric(particles.var.(sparse_fields{n_f})),'Sparse format currently only supports numeric variables');
                assert(size(particles.var.(sparse_fields{n_f}),2)==1,'In sparse format currently do not support variables being more than one dimensional');                
                particles.sparse_history.(sparse_fields{n_f})(:,n) = sparse(heritage,1,particles.var.(sparse_fields{n_f})(heritage),n_particles,1);
                particles.var.(sparse_fields{n_f}) = particles.var.(sparse_fields{n_f})(i_resample,:);
            end
            parents(:,n) = i_resample;
            for n_a = n:-1:2
                 if n_a~=1                    
                    n_times_sampled_local = accumarray(heritage, n_times_sampled_no_zeros, [n_particles,1], [], [], 1);
                    i_sampled = find(n_times_sampled_local);
                    n_times_sampled_no_zeros = nonzeros(n_times_sampled_local);
                    heritage = full(parents(i_sampled,n_a-1));
                 end
                
                i_sampled_reduced = parents(i_sampled,n_a-1);
                parents(:,n_a-1) = sparse(i_sampled,1,i_sampled_reduced,n_particles,1);
                for n_f = 1:numel(sparse_fields)
                    i_sampled_reduced = sort(i_sampled_reduced);
                    b_diff = [true;diff(i_sampled_reduced)~=0];
                    i_sampled_reduced = i_sampled_reduced(b_diff);
                    particles.sparse_history.(sparse_fields{n_f})(:,n_a-1) = ...
                        sparse(i_sampled_reduced,1,particles.sparse_history.(sparse_fields{n_f})(i_sampled_reduced,n_a-1),n_particles,1);
                end
            end
        else
            p_fields = fields(particles.var);
            for n_f = 1:numel(p_fields)
                particles.var.(p_fields{n_f}) = particles.var.(p_fields{n_f})(i_resample,:);
            end
        end
        
    end
end

z_max = max(log_weights);
w = exp(log_weights-z_max);
log_Z = log_Z+z_max+log(sum(w))-log(numel(w));
relative_weights = w/sum(w);

if b_sparse
    parents(:,end) = sparse((1:n_particles)');
    heritage = full(parents(:,end));    
    weights_no_zeros = relative_weights;
    weights_no_zeros(weights_no_zeros==0) = min(1e-32,min(weights_no_zeros(weights_no_zeros~=0))/1e10);
    sparse_relative_weights = sparse(n_particles,numel(sampling_functions));
    for n=numel(sampling_functions):-1:1
        weights_local = accumarray(heritage, weights_no_zeros, [n_particles,1], [], [], 1);
        i_sampled = find(weights_local);
        weights_no_zeros = nonzeros(weights_local);
        sparse_relative_weights(:,n) = accumarray(i_sampled, weights_no_zeros, [n_particles,1], [], [], 1);
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
    [particles, relative_weights] = compress_samples(particles,relative_weights,numel(weighting_functions));
end

if exist('sparse_relative_weights','var')
     particles.sparse_variable_relative_weights = sparse_relative_weights;
end

end
