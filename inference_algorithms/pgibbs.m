function [samples, log_Zs, mus, retained_particle] = pgibbs(sampling_functions,weighting_functions,...
               N,resample_method,n_iter,b_compress,f,b_Rao_Black,initial_retained_particle)
%pgibbs  Particle gibbs algorithm without static parameters
%
% Performs inference using iterated conditional SMC as described by section
% 2.2 in the paper.
%
% Required inputs:
%   sampling_functions = See infer.m
%   weighting_functions = See infer.m
%   N (+ve integer) = Number of particles, also N in paper
%   resample_method = Method used in resampling.  See resample_particles.m.
%                     If empty takes default from resample_particles.m
%   n_iter (+ve integer) = Number of iterations
%   b_compress (boolean) = Whether to use compress_samples
%   b_Rao_Black (boolean) = Whether to Rao-Blackwellize and return all
%                           generated samples or just the retained particle
%
% Optional inputs:
%   initial_retained_particle (stack_object) = Allows the algorithm to be
%                          initialized with a retained particle.  If not
%                          provided, the first iteration runs as an
%                          unconditional sweep.
%
% Outputs:
%   samples = Object array of type stack_object containing details about
%             sampled variables, their weights and any constant variables
%   log_Zs = Marginal likelihood of individual sweeps
%
% Tom Rainforth 07/06/16

log_Zs = NaN(n_iter,1);

if ~exist('initial_retained_particle','var'); initial_retained_particle = []; end
if ~exist('b_compress','var') || isempty(b_compress); b_compress = false; end
if ~exist('f','var'); f = []; end
if ~exist('b_Rao_Black','var') || isempty(b_Rao_Black); b_Rao_Black = true; end

b_compress = b_compress && b_Rao_Black;

retained_particle = initial_retained_particle;

for iter=1:n_iter
    [samples(iter), log_Zs(iter), retained_particle, mus(iter,:)] ... %FIXME preallocate
        = pg_sweep(sampling_functions,weighting_functions,N,retained_particle,resample_method,b_compress,f,b_Rao_Black); %#ok<AGROW>
        
    if iter==1 && ~b_compress
        % Memory management once have information from the first iteration
        [samples,b_compress] = memory_check(samples,n_iter,numel(sampling_functions));
        samples = repmat(samples,n_iter,1);
    end
    
end


end