function [samples, log_Zs, node_weights, sampled_indices, switching_rate, mus] = ...
    ipmcmc(sampling_functions,weighting_functions,N,resample_method,n_iter,b_compress,f,...
    b_Rao_Black,b_parallel,M,P,n_conditional_gibbs_cycles,initial_retained_particles)
%ipmcmc  iPMCMC inference algorithm
%
% Performs inference using interacting particle Markov chain Monte Carlo.
% For details please see the paper:
%   T. Rainforth, C. A. Naesseth, F. Lindsten, B. Paige, J-W. van de Meent,
%   A. Doucet and F. Wood. Interacting particle Markov chain Monte Carlo.
%   In ICML 2016
% Though can be used directly, it is recommended to use this through the 
% infer.m function which deals with many processes common to the different 
% inference algorithms such as allowing for default options and processing 
% the outputs.  In particular, infer converts model files to the
% automatically produce the sampling_functions and weighting_functions.
%
% Inputs:
%   sampling_functions = See infer.m
%   weighting_functions = See infer.m
%   N (+ve integer) = Number of particles, also N in paper
%   resample_method = Method used in resampling.  See resample_particles.m.
%                     If empty takes default from resample_particles.m
%   n_iter (+ve integer) = Number of MCMC iterations (indexed by r in the paper)
%   b_compress (boolean) = Whether to use compress_samples
%   f = Function to take expectation of.  Takes the var field of samples as
%       inputs.  See function_expectation.m.
%                               Default = []; (i.e. no estimate made)
%   b_Rao_Black (boolean | 'cond_update_only') = Controls level of
%                     Rao-Blackwellization (RB).  If true both the sampling of
%                     the retained particles and the conditional node
%                     updates are RBed, if 'cond_update_only' then only the
%                     latter is.  RB the sampling of the retained particle
%                     requires significantly more memory but should give
%                     more accurate results.  There is little loss in RB
%                     the conditional node update.
%   b_parrallel (boolean) = Whether to carry out the calculations in
%                           parallel
%   M (+ve integer) = Total number of nodes, also M in paper
%   P (+ve integer) = Number of conditional nodes
%   n_conditional_gibbs_cycles (+ve integer) = Number of gibbs sweeps to
%                     update the indices of the conditional nodes.  The
%                     main advantage of running more sweeps is to improve
%                     the estimate of the node weights for
%                     Rao-Blackwellization.  Typically 1 is sufficient.
%                           Default = 1
%   initial_retained_particles (stack_object array) = Allows the algorithm 
%                     to be initialized with a retained particle set.  If not
%                     provided, the first iteration for all nodes runs as an
%                     unconditional sweep.
%
% Outputs:
%   samples = Object array of type stack_object containing details about
%             sampled variables, their weights and any constant variables
%   log_Zs = Marginal likelihood of individual sweeps
%   node_weights = Relative weights given to each node under the RB at each
%            iterationb
%   sampled_indices = Indices of the conditional nodes for the NEXT
%            set of SMC sweeps (c_{1:P} in the paper).
%   switching_rate = Proportion of retained particles inhereted from an
%             unconditional sweep at each iteration
%   mus = Mean estimates of individual sweeps
%
% Tom Rainforth 07/06/16

if ~exist('f','var'); f = []; end

if ~exist('n_conditional_gibbs_cycles','var')
    n_conditional_gibbs_cycles = 1;
end

% Assign initial retained particles.  Whenever the retained particle is
% empty, pg_sweep will run an smc sweep rather than a csmc sweep.  Thus the
% cells of retained_particles corresponding to unconditional nodes are
% simply left empty.
if exist('initial_retained_particles','var') && ~isempty(initial_retained_particles)
    assert(P==numel(initial_retained_particles),'If providing initial_retained_particles, must provide P of them');
    retained_particles = [reshape(num2cell(initial_retained_particles),P,1);
        cell(M-P,1)];
else
    retained_particles = cell(M,1);
end

% Initialize the conditional node ids to the first P nodes
c = 1:P;
c_old = c;

% Preallocate space to some of stored outputs
[log_Zs,node_weights] = deal(NaN(n_iter,M));
sampled_indices = NaN(n_iter,P);
switching_rate = NaN(n_iter,1);
mus = cell(n_iter,M);

b_Rao_Black_particle_choice = ~strcmpi(b_Rao_Black,'cond_update_only') && b_Rao_Black;
b_Rao_Black_conditional_node_choice = strcmpi(b_Rao_Black,'cond_update_only') || b_Rao_Black;

b_compress = b_compress && b_Rao_Black_particle_choice;

for iter=1:n_iter
    
    % Local space preallocation required for parfor
    [particles_iter,new_retained_particle] = deal(cell(M,1));
    log_Zs_iter = NaN(M,1);
    mus_iter = cell(M,1);
    
    % Call the pg_sweep for each node.  When their is a retained particle
    % this is a csmc sweep, when there is none it is an ordinary
    % unconditional smc sweep.
    if b_parallel
        parfor node = 1:M
            [particles_iter{node},log_Zs_iter(node),new_retained_particle{node},mus_iter{node}] = ...
                pg_sweep(sampling_functions,weighting_functions,N,retained_particles{node},resample_method,false,f,b_Rao_Black_particle_choice);
        end
    else
        for node = 1:M
            [particles_iter{node},log_Zs_iter(node),new_retained_particle{node},mus_iter{node}] = ...
                pg_sweep(sampling_functions,weighting_functions,N,retained_particles{node},resample_method,false,f,b_Rao_Black_particle_choice);
        end
    end
    
    % Update the indices of the conditional nodes and establish the node
    % weightings required for Rao-Blackwellization
    
    % * Uncomment below and ** for checking against old code with fixed seed
    % c = 1:P;
    
    [c, node_weights_iter] = gibbs_update_c(log_Zs_iter(:),c,n_conditional_gibbs_cycles);
    retained_particles = cell(M,1);
    for m=1:P
        retained_particles{c(m)} = new_retained_particle{c(m)};
        % ** Uncomment below and * for checking against old code with fixed seed
        % retained_particles{m} = new_retained_particle{c(m)};
    end
        
    % Setup output samples
    for m=1:M
        if b_Rao_Black_conditional_node_choice
            % Rao Blackwellize the Gibbs update of the conditional nodes
            samples(iter,m) = particles_iter{m}; %#ok<AGROW>
            samples(iter,m).relative_particle_weights = samples(iter,m).relative_particle_weights.*node_weights_iter(m); %#ok<AGROW>
            if b_compress && b_Rao_Black_particle_choice
                samples(iter,m) = compress_samples(samples(iter,m),numel(weighting_functions)); %#ok<AGROW>
            end
        elseif m<=P
            % No Rao Blackwellization. Add when m<=P.
            samples(iter,m) = retained_particles{c(m)}; %#ok<AGROW>
            samples(iter,m).relative_particle_weights = 1/P; %#ok<AGROW>
        else
            continue
        end
    end
    
    % Store additional outputs outputs
    log_Zs(iter,:) = log_Zs_iter;
    mus(iter,:) = mus_iter;
    node_weights(iter,:) = node_weights_iter;
    sampled_indices(iter,:) = c;
    switching_rate(iter) = (P-numel(intersect(c,c_old)))/P;
    c_old = c;
    
    if iter==1 && ~b_compress
        % Memory management once have information from the first iteration
        [samples,b_compress] = memory_check(samples,n_iter,numel(sampling_functions));
        % Preallocate space for the other samples
        samples = repmat(samples,n_iter,1);
    end
       
end



end