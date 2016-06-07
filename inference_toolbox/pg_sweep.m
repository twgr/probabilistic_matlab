function [particles, log_Z, retained_particle] = pg_sweep(sampling_functions,...
                weighting_functions,N,retained_particle,b_compress,b_Rao_Black)
%pg_sweep   Sweep used for Particle Gibbs and iPMCMC
%
% When provided with a retained particle, performs a conditional sequential
% Monte Carlo (CSMC) sweep as per Algorithm 2 in the paper.  Otherwise
% Performs sequential Monte Carlo (SMC) as per Algorithm 1 in the paper.
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
%   b_compress (boolean) = Exploits the degeneracy caused by resampling to
%                          store the output using sparse matrices and an
%                          expliticly stored, sparse, ancestral lineage.
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
%
% Tom Rainforth 07/06/16

global sample_size

if ~exist('retained_particle','var') || isempty(retained_particle) || isnan(retained_particle)
    % Run as an unconditional smc sweep if no retained particle supplied
    [particles, log_Z] = smc_sweep(sampling_functions,weighting_functions,N,false,false);
else
    % Otherwise run as a conditional smc sweep.
    
    % Global variable that can be used as a controller inside the code for
    % sampling_functions and weighting functions if desired    
    sample_size = N-1; %#ok<NASGU>    
    particles = stack_object;
    log_Z = 0;
    
    for n=1:numel(sampling_functions)
        % Sample from eq 1a in the paper
        particles = sampling_functions{n}(particles);        
        log_weights = [weighting_functions{n}(retained_particle); % Retained particle weight
                       weighting_functions{n}(particles)];        % Other particle weights
        
        if n~=numel(sampling_functions)
            % At each step except the last, perform the conditional
            % resampling step of the particles.  This is broken up into
            % resampling of all the variables rather than treating seperate
            % particles as fully seperate indices for speed purposes.
            
            % Resample indices
            [i_resample, log_Z_step] = resample_step(log_weights, numel(log_weights)-1);
            log_Z = log_Z+log_Z_step;
            
            % Particle fields that are being resampled 
            p_fields = fields(particles.var);
            for n_f = 1:numel(p_fields)
                % For each field, the sampled particles are the union of
                % the resampled particles originating for the unretained
                % particles and sum number of replications of the retained
                % particle.  Note that we are generating N-1 samples here,
                % with the retained particle completing the set
                particles.var.(p_fields{n_f}) = [particles.var.(p_fields{n_f})(i_resample(i_resample~=1)-1,:);...
                                                 repmat(retained_particle.var.(p_fields{n_f})(1,1:size(particles.var.(p_fields{n_f}),2)),sum(i_resample==1),1)];
            end
        end
    end
    
    % Add the retained particle back into the final particle set
    p_fields = fields(particles.var);
    for n_f = 1:numel(p_fields)
        particles.var.(p_fields{n_f}) = [retained_particle.var.(p_fields{n_f});...
                                         particles.var.(p_fields{n_f})];
    end
    
    % Calculate the marginal likelihood and the relative particle weights
    z_max = max(log_weights);
    w = exp(log_weights-z_max);
    log_Z = log_Z+z_max+log(sum(w))-log(numel(w));
    particles.relative_particle_weights = w/sum(w);    
end

% Sample the new retained particle
i_keep = datasample(1:numel(particles.relative_particle_weights),1,'Weights',particles.relative_particle_weights,'Replace',true);
retained_particle = stack_object;
retained_particle.con = particles.con;
for n_f = 1:numel(p_fields)
    retained_particle.var.(p_fields{n_f}) = particles.var.(p_fields{n_f})(i_keep,:);
end

% Reset the sample_size variable
sample_size = N;

if ~b_Rao_Black
    particles = retained_particle;
    particles.relative_particle_weights = 1;
elseif b_compress
    particles = compress_samples(particles,numel(weighting_functions));
end

end