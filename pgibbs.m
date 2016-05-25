function [samples, relative_weights, log_Zs] = pgibbs(sampling_functions,weighting_functions,n_particles,n_iter,b_compress,initialization)

relative_weights = NaN(n_particles*n_iter,1);
log_Zs = NaN(n_iter,1);

if isempty(initialization)
    [samples, log_Zs(1), relative_weights(1:n_particles)] = smc_sweep(sampling_functions,weighting_functions,n_particles,false); % FIXME add ability to use sparse to pgiibs
    retained_particle = stack_object;
    retained_particle.con = samples.con;
    p_fields = fields(samples.var);
    i_keep = datasample(1:n_particles,1,'Weights',relative_weights(1:n_particles,:),'Replace',true);
    
    for n_f = 1:numel(p_fields)
        retained_particle.var.(p_fields{n_f}) = samples.var.(p_fields{n_f})(i_keep,:);
    end
else
    [samples, log_Zs(1), relative_weights(1:n_particles), retained_particle] ...
        = pg_sweep(sampling_functions,weighting_functions,n_particles,initialization,false);
end

if b_compress
    [samples,relative_weights(1:n_particles)] = compress_samples(samples,relative_weights(1:n_particles),numel(weighting_functions));
end

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
            [samples, relative_weights(1:n_particles)] = compress_samples(samples, relative_weights(1:n_particles), numel(sampling_functions));
        end
    end
end

samples = repmat(samples,n_iter,1);

for n=2:n_iter
    [samples(n), log_Zs(n), ...
        relative_weights((1+(n-1)*n_particles):(n*n_particles)), retained_particle] ...
        = pg_sweep(sampling_functions,weighting_functions,n_particles,retained_particle,b_compress);
end


end