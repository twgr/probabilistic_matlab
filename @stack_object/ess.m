function E = ess(samples,field,dims,i_samples)
%ess Effective sample size calculation
%
% E = ess(samples,field,dims,i_samples)
%
% Calculates the effective sample size as defined in the ipmcmc paper.
%
% Inputs:
%   samples = stack_object to calculate ESS for
%   field = String specifying required variable field (only one accepted at
%           a time).
%   dims = Dimensions to calculate the ess for.  By default, all are
%          calculated.
%   i_samples = The indices of the samples to use in the calculate.  By
%               default, all are used.
%
% Outputs:
%   E = Effective sample size.  Returned as 1xd vector giving the ESS for
%       each dimension.
%
% Tom Rainforth 08/07/16

[n_samples,n_dims_total] = size(samples.var.(field));

if ~exist('dims','var') || isempty(dims)
    dims = 1:n_dims_total;
end

if ~exist('i_samples','var') || isempty(i_samples)
    i_samples = (1:n_samples)';
end

E = NaN(1,numel(dims));

for nd = 1:numel(dims)
    d=dims(nd);
    x_local = samples.get_variable_dim(field,d,i_samples);
    w_local = samples.get_weights(field,d,i_samples);
    [~,~,i_val] = fast_unique(x_local);
    V = accumarray(i_val,w_local);
    E(nd) = 1./sum(V.^2);
end

end