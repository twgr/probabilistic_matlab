function [samples, log_Zs, b_accept] = a_pgibbs(sampling_functions,weighting_functions,...
    N,n_iter,b_compress,b_Rao_Black,initial_retained_particle)
%a_pgibbs   Alternate move particle Gibbs 
%
% Carries out the alternate move particle Gibbs (APG) algorithm which
% interleaves PG and PIMH steps.  For more information see section 4 of the
% iPMCMC paper or Roman Holenstein's PhD thesis.
%
% Required inputs:
%   sampling_functions = See infer.m
%   weighting_functions = See infer.m
%   N (+ve integer) = Number of particles, also N in paper
%   n_iter (+ve integer) = Number of iterations
%   b_Rao_Black (boolean) = Whether to Rao-Blackwellize and return all
%                           generated samples or just the retained particle
%   b_compress (boolean) = Whether to use compress_samples
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
%   b_accept = Boolean vector indicating if that iteration is accepted
%
% Tom Rainforth 08/06/16


log_Zs = NaN(n_iter,1);
b_accept = true(n_iter,1);
b_compress = b_compress && b_Rao_Black;

if ~exist('initial_retained_particle','var')
    initial_retained_particle = [];
end
retained_particle = initial_retained_particle;

for iter=1:n_iter
    
    if mod(iter,2)==0
        % At even steps, carry out a PIMH transition.  This uses pg_sweep
        % rather than smc_sweep because of requiring the sampling of a new
        % retained particle at the end of the sweep
        [sample_proposed, log_Z_proposed, retained_particle_proposed] ...
            = pg_sweep(sampling_functions,weighting_functions,N,[],b_compress,b_Rao_Black);
        bKeep = rand<min(1,exp(log_Z_proposed-log_Zs(iter-1)));
        if bKeep
            samples(iter) = sample_proposed; %#ok<AGROW>
            log_Zs(iter) = log_Z_proposed;
            retained_particle = retained_particle_proposed;
        else
            samples(iter) = samples(iter-1); %#ok<AGROW>
            log_Zs(iter) = log_Zs(iter-1);
            b_accept(iter) = false;
        end
    else
        % At the odd steps, carry out a PG step
        [samples(iter), log_Zs(iter), retained_particle] ...
            = pg_sweep(sampling_functions,weighting_functions,N,retained_particle,b_compress,b_Rao_Black); %#ok<AGROW>
    end
    
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