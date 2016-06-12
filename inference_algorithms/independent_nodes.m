function [samples, log_Zs, b_accepts] = independent_nodes(sampling_functions,weighting_functions,...
                N,n_iter,b_compress,b_Rao_Black,b_parallel,Ms,initial_retained_particles)
%independent_nodes  Multiple non-interacting PMCMC nodes
%
% Allows an arbitrary combination of non-interacting particle gibbs, PIMH
% and alternative move particle gibbs nodes.  Can be used to generate
% results for mPG, mPIMH and mAPG as per the paper in additition to other
% combinations (e.g. the half-half results in the supplementary material
% using M/2 PG nodes and M/2 PIMH nodes).
%
% Inputs:
%   sampling_functions = See infer.m
%   weighting_functions = See infer.m
%   N (+ve integer) = Number of particles, also N in paper
%   n_iter = Number of independent sweeps to perform
%   b_compress (boolean) = Whether to use compress_samples
%   b_Rao_Black (boolean) = Whether to Rao-Blackwellize and return all
%                           generated samples or just the retained particle
%   b_parallel (boolean) = Whether to run sweeps in parallel (number of
%                          threads automatically set by parfor)
%   Ms (3x1 vector of +ve integers) = Number of nodes dedicated to
%                          particle Gibbs, alternate move PG and PIMH
%                          respectively.  For example, [0,32,0] corresponds
%                          to mPIMH in the paper with M=32.
%   initial_retained_particles (stack_object array) = Allows the algorithm 
%                     to be initialized with a retained particle set.  If not
%                     provided, the first iteration for all nodes runs as an
%                     unconditional sweep.
%
% Output:
%   samples = Object array of type stack_object containing details about
%             sampled variables, their weights and any constant variables
%   log_Zs = Marginal likelihood of individual sweeps
%   b_accept = Boolean vector indicating if that iteration is accepted for
%            PIMH / PIMH steps of M_apg nodes.  NaN indicates that the
%            corresponding step / node is not a PIMH iteration.
%
% Tom Rainforth 12/06/16

M_pg = Ms(1);
M_apg = Ms(2);
M_pimh = Ms(3);

M = M_pg+M_apg+M_pimh;

samples = repmat(stack_object,n_iter,M);
log_Zs = NaN(n_iter,M);
b_accepts = NaN(n_iter,M);

% Note that parfor needs an instance in each branch, even though most are
% irrelevant
if isempty(initial_retained_particles)
    initial_retained_particles = cell(M,1);
elseif ~iscell(initial_retained_particles)
    initial_retained_particles = num2cell(initial_retained_particles(:),2);
    initial_retained_particles = [initial_retained_particles; repmat(initial_retained_particles,M-size(initial_retained_particles,1),1)];
end

if b_parallel    
    parfor n_n = 1:M
        if n_n<=M_pg
            [samples(:,n_n), log_Zs(:,n_n)] = pgibbs(sampling_functions,weighting_functions,N,n_iter,b_compress,b_Rao_Black,initial_retained_particles{n_n});
        elseif n_n<=(M_pg+M_apg)
            [samples(:,n_n), log_Zs(:,n_n),b_accepts(:,n_n)] = a_pgibbs(sampling_functions,weighting_functions,N,n_iter,b_compress,b_Rao_Black,initial_retained_particles{n_n});
        else
            [samples(:,n_n), log_Zs(:,n_n),b_accepts(:,n_n)] = pimh(sampling_functions,weighting_functions,N,n_iter,b_compress,b_Rao_Black);
        end
    end    
else    
    for n_n = 1:M
        if n_n<=M_pg
            [samples(:,n_n), log_Zs(:,n_n)] = pgibbs(sampling_functions,weighting_functions,N,n_iter,b_compress,initial_retained_particles{n_n});
        elseif n_n<=(M_pg+M_apg)
            [samples(:,n_n), log_Zs(:,n_n),b_accepts(:,n_n)] = a_pgibbs(sampling_functions,weighting_functions,N,n_iter,b_compress,b_Rao_Black,initial_retained_particles{n_n});
        else
            [samples(:,n_n), log_Zs(:,n_n),b_accepts(:,n_n)] = pimh(sampling_functions,weighting_functions,N,n_iter,b_compress,b_Rao_Black);
        end
    end    
end

end



