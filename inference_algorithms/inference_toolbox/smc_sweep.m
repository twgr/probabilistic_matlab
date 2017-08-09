function [particles, log_Z, mu] = smc_sweep(sampling_functions,...
                        weighting_functions,N,resample_method,b_compress,f,prop_sub_sample)
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
%   resample_method = Method used in resampling.  See resample_particles.m.
%                     If empty takes default from resample_particles.m
%                               Default = []
%   b_compress (boolean) = Whether to use compress_samples
%                               Default = false;
%   f = Function to take expectation of.  Takes the var field of samples as
%       inputs.  See function_expectation.m.
%                               Default = []; (i.e. no estimate made)
%   prop_sub_sample (0<x<=1) = Proportion of samples to return.  For memory
%                     reasons it may be beneficial to subsample the 
%                     produced particles down to a smaller number (similar
%                     to not Rao-Blackwellizing for pmcmc methods).
%                               Default = 1;
%
% Outputs:
%   particles = Object of type stack_object storing all the samples.
%   log_Z = Log marginal likelihood estimate as per equation 4 in the paper
%   mu = Expectation for this sweep
%
% Tom Rainforth 07/06/16

if ~exist('resample_method','var'); resample_method = []; end
if ~exist('b_compress','var') || isempty(b_compress); b_compress = false; end
if ~exist('prop_sub_sample','var') || isempty(prop_sub_sample); prop_sub_sample = 1; end
if ~exist('f','var'); f = []; end

% Global parameter that allows for communication of the sample size with
% the sampling_functions and weighting_functions if desired.
global sample_size
sample_size = N;

particles = stack_object;
log_Z = 0;

for n=1:numel(sampling_functions)
    % Simulate forwards as per eq 1a in the paper
    particles = sampling_functions{n}(particles);
    % Weight as per eq 1b in the paper
    log_weights = weighting_functions{n}(particles);
    
    if n~=numel(sampling_functions)
        % When not at the final step, perform resampling.
        [particles, log_Z_step] = resample_particles(particles, log_weights, N, resample_method);
        log_Z = log_Z+log_Z_step;
    end
end

% Set outputs
z_max = max(log_weights);
w = exp(log_weights-z_max);
log_Z = log_Z+z_max+log(sum(w))-log(numel(w)); % Marginal likelihood
particles.relative_particle_weights = w/sum(w);

if prop_sub_sample~=1
    % Sub sample down as requeseted
     particles = resample_particles(particles, log_weights, N*prop_sub_sample, resample_method);
end

% Calculate function expectations
if ~isempty(f)
    mu = particles.function_expectation(f);
else
    mu = zeros(1,0);
end

% Tying up of sample compression as required.
if b_compress
    particles = compress_samples(particles,numel(weighting_functions));
end


end
