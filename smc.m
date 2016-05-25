function [samples, relative_weights, log_Z_total, log_Zs] = smc(sampling_functions,weighting_functions,n_particles,n_iter,b_compress,b_sparse)

if b_sparse && b_compress
    warning('b_sparse does its own compression, no extra compression required');
    b_compress = false;
end

log_Zs = NaN(n_iter,1);
relative_weights_n = cell(n_iter,1);
[samples, log_Zs(1), relative_weights_n{1}] = smc_sweep(sampling_functions,weighting_functions,n_particles,b_compress,b_sparse);

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
            [samples, relative_weights_n{1}] = compress_samples(samples, relative_weights_n{1}, numel(sampling_functions));
        end
    end
end

samples = repmat(samples,n_iter,1);

for n=2:n_iter
    [samples(n), log_Zs(n), relative_weights_n{n}] = smc_sweep(sampling_functions,weighting_functions,n_particles,b_compress,b_sparse);
end

max_log_Z = max(log_Zs);
w = exp(log_Zs-max_log_Z);

log_Z_total = log(mean(w))+max_log_Z;

w = w/sum(w);
relative_weights = NaN(n_particles,1);
for n=1:n_iter
    relative_weights((1+(n-1)*n_particles):(n*n_particles)) = w(n)*relative_weights_n{n};
end

if ~isempty(samples(1).sparse_variable_relative_weights)
    if isstruct(samples(1).sparse_variable_relative_weights)
        sparse_fields = fields(samples(1).sparse_variable_relative_weights);
        for n=1:n_iter
            for n_f=1:numel(sparse_fields)
                samples(n).sparse_variable_relative_weights.(sparse_fields{n_f}) = samples(n).sparse_variable_relative_weights.(sparse_fields{n_f})*w(n);
            end
        end
    else
        for n=1:n_iter
            samples(n).sparse_variable_relative_weights = samples(n).sparse_variable_relative_weights*w(n);
        end
    end
end

end