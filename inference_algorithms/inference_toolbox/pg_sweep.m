function [particles, log_Z, retained_particle, mu] = pg_sweep(sampling_functions,...
    weighting_functions,N,retained_particle,resample_method,b_compress,f,b_Rao_Black)
%pg_sweep   Sweep used for Particle Gibbs, alternative move PG and iPMCMC
%
% When provided with a retained particle, performs a conditional sequential
% Monte Carlo (CSMC) sweep as per Algorithm 2 in the paper.  Otherwise
% Performs sequential Monte Carlo (SMC) as per Algorithm 1 in the paper.
% In the latter case, behaviour differs from smc_sweep by sampling a
% retained particle at the end of the sweep and storing required
% information to reassemble intermediary retained particles.
%
% Required inputs:
%   sampling_functions = See infer.m
%   weighting_functions = See infer.m
%   N (+ve integer) = Number of particles, also N in paper
%
% Optional inputs:
%   retained_particle = A retained particle of type stack_object.  Sweep
%                       will be conditioned on this particle when provided.
%                               Default = empty
%   resample_method = Method used in resampling.  See resample_particles.m.
%                     If empty takes default from resample_particles
%                               Default = []
%   b_compress (boolean) = Whether to use compress_samples
%                               Default = false;
%   b_Rao_Black (boolean) = If true all particles are returned with
%                           weights, else only a single particle (i.e. the
%                           retained particle) is returned
%                               Default = true;
%
% Outputs:
%   particles = Object of type stack_object storing all the samples.
%   log_Z = Log marginal likelihood estimate as per equation 4 in the paper
%   retained_particle = Particle sampled as the retained particle as per
%                       equation 7 (sampled in proportion to weight)
%   mu = Expectation for this sweep
%
% Tom Rainforth 07/06/16

global sample_size

% Set defaults on anything not provided
if ~exist('retained_particle','var'); retained_particle = []; end
if ~exist('resample_method','var'); resample_method = []; end
if ~exist('b_compress','var') || isempty(b_compress); b_compress = false; end
if ~exist('f','var'); f = []; end
if ~exist('b_Rao_Black','var') || isempty(b_Rao_Black); b_Rao_Black = true; end

b_conditional_sweep = ~isempty(retained_particle);

% Global variable that can be used as a controller inside the code for
% sampling_functions and weighting functions if desired
if b_conditional_sweep
    sample_size = N-1; %#ok<NASGU>
    n_resample = N-1;
else
    sample_size = N; %#ok<NASGU>
    n_resample = N;
end

particles = stack_object;
log_Z = 0;
% These variables store intermediary details allows for later
% reconstruction
[variables_step,sizes_step] = deal(cell(numel(sampling_functions),1));

for n=1:numel(sampling_functions)
    % Sample from eq 1a in the paper
    particles = sampling_functions{n}(particles);
    
    % Store the variables at this step and the sizes for allowing later
    % reconstruction of the intermediary retained particles from the
    % final retained particle.
    variables_step{n} = fields(particles.var);
    for v=1:numel(variables_step{n})
        sizes_step{n}{v} = size(particles.var.(variables_step{n}{v}),2);
    end
    
    if b_conditional_sweep
        % For the conditional then we need to add the retained particle
        % back into the particle set.  This is done by reconstructing the
        % intermediary retained particle using the previously stored
        % information
        
        % Reconstruction intermediary retained particle
        intermediary_retained_particle = stack_object;
        intermediary_retained_particle.con = particles.con;
        for v=1:numel(retained_particle.variables_step{n})
            if ~iscell(retained_particle.var.(retained_particle.variables_step{n}{v}))
                intermediary_retained_particle.var.(retained_particle.variables_step{n}{v}) = ...
                    retained_particle.var.(retained_particle.variables_step{n}{v})(1:retained_particle.sizes_step{n}{v});
            else
                intermediary_retained_particle.var.(retained_particle.variables_step{n}{v}) = ...
                    {retained_particle.var.(retained_particle.variables_step{n}{v}){1}(1:retained_particle.sizes_step{n}{v})};
            end
        end
                
        % Add the intermediary retained particle back into the particle set
        % TODO - the below is an elegant but inefficient way of doing the
        % conditional resampling and currently adds ~20% to the run time.
        % Switch to using a conditional resample directly except at the
        % final iteration.
        particles = compose_two_sample_objects(intermediary_retained_particle,particles,':',1,N-1);
        
    end
    
    % Calculate the log weight
    log_weights = weighting_functions{n}(particles);
    
    if n~=numel(sampling_functions)
        % Resample indices.  Note for conditional that this does not include the
        % retained particle that is forced to survive, but can include
        % additional instances generated by the resampling.
        [particles,log_Z_step] = resample_particles(particles, log_weights, n_resample, resample_method);
        log_Z = log_Z+log_Z_step;
    end
end

% Sample the new retained particle (this is equivalent to resampling down
% to a single particle)
[retained_particle, log_Z_step] = resample_particles(particles, log_weights, 1);

% Update log_Z estimate for last step
log_Z = log_Z+log_Z_step;

% Store the variables that exist at each step and their sizes
retained_particle.variables_step = variables_step;
retained_particle.sizes_step = sizes_step;

z_max = max(log_weights);
w = exp(log_weights-z_max);
particles.relative_particle_weights = w/sum(w);

% Reset the sample_size variable
sample_size = N;

% FIXME, might not want to calculate f at all points as could be expensive,
% should add the option to calculate f after the rao blackwellisation
if ~isempty(f)
    mu = particles.function_expectation(f);
else
    mu = zeros(1,0);
end

% If not Rao Blackwellizing, the returned particle is the retained
% particle rather than the full set.
if ~b_Rao_Black
    particles = retained_particle;
    particles.relative_particle_weights = 1;
end

if b_compress
    particles = compress_samples(particles,numel(weighting_functions));
end

end