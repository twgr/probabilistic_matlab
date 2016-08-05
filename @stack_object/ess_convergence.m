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

E = NaN(n_iter,numel(dims));

for nd = 1:numel(dims)
    d=dims(nd);
    x_local = samples.get_variable(field,d);
    [w_local,~,iNonZeros] = samples.get_weights(field,d);
    i_iter = ceil(iNonZeros/(n_p_per_iter));  
    [~,~,i_val] = fast_unique(x_local);
    V = accumarray([i_iter,i_val],w_local,[n_iter,max(i_val)],[],[],1)';
    % Matlab is quicker working with columns of sparse arrays than rows
    n_vi = size(V,1);
    vi = zeros(n_vi,1);
    for ni = 1:size(V,2)
        ithis = find(V(:,ni));
        vi_old_this = vi(ithis);
        vi(ithis) = vi_old_this+nonzeros(V(:,ni));
        % Would ideally not want to do this in the loop, but full(V) may be
        % a larger matrix than can be stored
        if ni==1 || (numel(ithis) > n_vi/2)
            a = sum(vi);
            a2 = sum(vi.^2);
        else
            a = (a+sum(nonzeros(V(:,ni))));
            a2 = (a2+sum(vi(ithis).^2)-sum(vi_old_this.^2));
        end
        E(ni,nd) = a^2/a2;
    end
end

end