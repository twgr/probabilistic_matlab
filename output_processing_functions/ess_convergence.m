function E = ess_convergence(samples,field,dims)
%ess Effective sample size convergence calculation
%
% E = ess_convergence(samples,field,dims)
%
% Calculates the effective sample size as defined in the ipmcmc paper after
% each mcmc iteration (or sweep in the case of smc
%
% Inputs:
%   samples = stack_object to calculate ESS for
%   field = String specifying required variable field (only one accepted at
%           a time).
%   dims = Dimensions to calculate the ess for.  By default, all are
%          calculated.
%
% Outputs:
%   E = Effective sample size.  Returned as n_iterxd vector giving the ESS 
%       for each dimension and iteration.
%
% Tom Rainforth 08/07/16

x = samples.var.(field);

if ~exist('dims','var')
    dims = 1:size(x,2);
end

n_iter = samples.options.n_iter;
n_p_per_iter = size(x,1)/n_iter;

if ~isempty(samples.relative_particle_weights)
    w_particles = samples.relative_particle_weights;
else
    w_particles = ones(numel(i_samples),1)/numel(i_samples);
end

E = NaN(n_iter,numel(dims));

for d=dims
    if issparse(x)
        iNonZeros = find(x(:,d));
        x_local = nonzeros(x(:,d));
    else
        iNonZeros = 1:size(x_local,2);
        x_local = x(:,d);
    end
    if isempty(samples.sparse_variable_relative_weights)
        % Includes support of outdated form
        w_local = w_particles;
    elseif isnumeric(samples.sparse_variable_relative_weights)
        w_local = nonzeros(samples.sparse_variable_relative_weights(:,d));
    else
        w_local = nonzeros(samples.sparse_variable_relative_weights.(field)(:,d));
    end
    
    i_iter = ceil(iNonZeros/(n_p_per_iter));  
    w_local = w_local/sum(w_local);
    [~,~,i_val] = fast_unique(x_local);
    V = accumarray([i_iter,i_val],w_local,[n_iter,max(i_val)],[],[],1);
    V = cumsum(V,1);
    E(:,d) = sum(V,2).^2./sum(V.^2,2);
end

end