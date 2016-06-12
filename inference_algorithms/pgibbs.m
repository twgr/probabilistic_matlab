function [samples, log_Zs] = pgibbs(sampling_functions,weighting_functions,...
                                    N,n_iter,b_compress,b_Rao_Black,initial_retained_particle)
%pgibbs  Particle gibbs algorithm without static parameters
%
% Performs inference using iterated conditional SMC as described by section
% 2.2 in the paper.
%
% Required inputs:
%   sampling_functions = See infer.m
%   weighting_functions = See infer.m
%   N (+ve integer) = Number of particles, also N in paper
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

if ~exist('initial_retained_particle','var')
    initial_retained_particle = [];
end

b_compress = b_compress && b_Rao_Black;

retained_particle = initial_retained_particle;

for iter=1:n_iter
    [samples(iter), log_Zs(iter), retained_particle] ...
        = pg_sweep(sampling_functions,weighting_functions,N,retained_particle,b_compress,b_Rao_Black); %#ok<AGROW>
        
    if iter==1
        %% Memory management once have information from the first iteration
        if b_Rao_Black && ~b_compress 
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
                    samples = compress_samples(samples, numel(sampling_functions));
                end
            end
        end
        
        samples = repmat(samples,n_iter,1);
    end
    
end


end