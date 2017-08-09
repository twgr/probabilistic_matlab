function [particles, log_Z_step, i_resample, n_times_sampled] = ...
                resample_particles(particles,log_weights,n_samples,method)
%resample_particles   Resampling step of particles
%
% Given particles, weights, a number of samples and a method, performs the
% particle resampling step
%
% Inputs:
%   particles (stack_object) = Particles to be resampled
%   log_weights (vector of +ve reals) = Need not be normalized
%   n_samples (+ve integer) = Number of samples to take
%   method ('stratified' | 'systematic' | 'multinomial' | 'residual')
%       = Employed resampling method.  If empty takes default from 
%         resample_indices
%
% Outputs
%   particles = Resampled particles
%   log_Z_step = Log marginal likelihood of step
%   i_resample = Indices of the particles to keep
%   n_times_sampled = Number of times each particle was sampled
%
% Tom Rainforth 13/06/16

if ~exist('method','var') || isempty(method)
    method = [];
end

assert(~particles.compressed,...
        'Cannot use the resample_particles function on an object that \n has already been compressed.  Use weighted_samples_to_unweighted instead / first');

% Guard for numerical instability
n_samples = round(n_samples);

z_max = max(log_weights);
w = exp(log_weights-z_max);
log_Z_step = z_max+log(sum(w))-log(numel(w));

% Establish the resampled indices
[i_resample,n_times_sampled] = resample_indices(w,n_samples,method);

% Update particles object as per the resampling indices
p_fields = fields(particles.var);
for n_f = 1:numel(p_fields)
    particles.var.(p_fields{n_f}) = particles.var.(p_fields{n_f})(i_resample,:);
end

if ~isempty(particles.relative_particle_weights)
    [~,rel_weights] = log_sum_exp(log(particles.relative_particle_weights(i_resample,:))-log_weights(i_resample,:));
    particles.relative_particle_weights = rel_weights./sum(rel_weights);
end

end
