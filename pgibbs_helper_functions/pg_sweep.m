function [particles, log_Z, relative_weights, retained_particle] = pg_sweep(sampling_functions,weighting_functions,n_particles,retained_particle,b_compress)

global sample_size
sample_size = n_particles-1;

particles = stack_object;
log_Z = 0;

for n=1:numel(sampling_functions)
    particles = sampling_functions{n}(particles);
    if isempty(retained_particle.con)
        retained_particle.con = particles.con;
    end
    log_weights = [weighting_functions{n}(retained_particle);weighting_functions{n}(particles)];
    if n~=numel(sampling_functions)
        [i_resample, log_Z_step] = resample_step(log_weights, numel(log_weights)-1);
        log_Z = log_Z+log_Z_step;
        
        p_fields = fields(particles.var);
        for n_f = 1:numel(p_fields)
            particles.var.(p_fields{n_f}) = [particles.var.(p_fields{n_f})(i_resample(i_resample~=1)-1,:);...
                                             repmat(retained_particle.var.(p_fields{n_f})(1,1:size(particles.var.(p_fields{n_f}),2)),sum(i_resample==1),1)];
        end        
    end
end

p_fields = fields(particles.var);
for n_f = 1:numel(p_fields)
    particles.var.(p_fields{n_f}) = [retained_particle.var.(p_fields{n_f});...
                                     particles.var.(p_fields{n_f})];
end

z_max = max(log_weights);
w = exp(log_weights-z_max);
log_Z = log_Z+z_max+log(sum(w))-log(numel(w));
relative_weights = w/sum(w);

i_keep = datasample(1:numel(relative_weights),1,'Weights',relative_weights,'Replace',true);

retained_particle = stack_object;
retained_particle.con = particles.con;

for n_f = 1:numel(p_fields)
    retained_particle.var.(p_fields{n_f}) = particles.var.(p_fields{n_f})(i_keep,:);
end

sample_size = n_particles;

if b_compress
    [particles, relative_weights] = compress_samples(particles,relative_weights,numel(weighting_functions));
end

end