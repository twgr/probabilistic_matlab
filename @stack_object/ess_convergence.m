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

[n_samples,n_dims_total] = size(samples.var.(field));

if ~exist('dims','var') || isempty(dims)
    dims = 1:n_dims_total;
end

n_iter = samples.options.n_iter;
n_p_per_iter = n_samples/n_iter;

i_samples = (1:n_samples)';

E = NaN(n_iter,numel(dims));

for nd = 1:numel(dims)
    d=dims(nd);
    [x_local,iNonZeros] = samples.get_variable_dim(field,d,i_samples);
    w_local = samples.get_weights(field,d,i_samples);
    i_iter = ceil(iNonZeros/(n_p_per_iter));  
    [~,~,i_val] = fast_unique(x_local);
    V = accumarray([i_iter,i_val],w_local,[n_iter,max(i_val)],[],[],1);
    V = cumsum(V,1);
    E(:,nd) = sum(V,2).^2./sum(V.^2,2);
end

end