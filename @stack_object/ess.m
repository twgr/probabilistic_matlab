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

x = samples.var.(field);
n_samples_total = size(x,1);

if ~exist('dims','var')
    dims = 1:size(x,2);
end

if ~exist('i_samples','var') || isempty(i_samples)
    i_samples = (1:n_samples_total)';
end

if ~isempty(samples.relative_particle_weights)
    w_particles = samples.relative_particle_weights(i_samples,:);
else
    w_particles = ones(numel(i_samples),1)/numel(i_samples);
end

x = x(i_samples,:);
if isnumeric(samples.sparse_variable_relative_weights) && ~isempty(samples.sparse_variable_relative_weights)
    samples.sparse_variable_relative_weights = samples.sparse_variable_relative_weights(i_samples,:);
elseif isstruct(samples.sparse_variable_relative_weights)
    samples.sparse_variable_relative_weights.(field) = samples.sparse_variable_relative_weights.(field)(i_samples,:);
end

E = NaN(1,numel(dims));

for d=dims
    if issparse(x)
        x_local = nonzeros(x(:,d));
    else
        x_local = x(:,d);
    end
    if isempty(samples.sparse_variable_relative_weights)
        w_local = w_particles;
    elseif isnumeric(samples.sparse_variable_relative_weights)
        w_local = nonzeros(samples.sparse_variable_relative_weights(:,d));
    else
        w_local = nonzeros(samples.sparse_variable_relative_weights.(field)(:,d));
    end
    w_local = w_local/sum(w_local);
    [~,~,i_val] = fast_unique(x_local);
    V = accumarray(i_val,w_local);
    E(d) = 1./sum(V.^2);
end

end