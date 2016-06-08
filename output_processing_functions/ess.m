function E = ess(samples,field,i_samples)

x = samples.var.(field);

n_samples_total = size(x,1);
n_dims = size(x,2);

if ~exist('i_samples','var') || isempty(i_samples)
    i_samples = (1:n_samples_total)';
end

if isfield(samples,'relative_weights') && ~isempty(samples.relative_weights')
    % Support of outdated form
    w_particles = samples.relative_weights(i_samples,:);
    w_particles = w_particles/sum(w_particles);
elseif isprop(samples,'relative_particle_weights')  && ~isempty(samples.relative_particle_weights)
    w_particles = samples.relative_particle_weights(i_samples,:);
else
    w_particles = ones(numel(i_samples),1)/numel(i_samples);
end

x = x(i_samples,:);
if isfield(samples.var,'num_instances')
    samples.var.num_instances = samples.var.num_instances(i_samples,:);
elseif isprop(samples,'sparse_variable_relative_weights') && isnumeric(samples.sparse_variable_relative_weights) && ~isempty(samples.sparse_variable_relative_weights)
    samples.sparse_variable_relative_weights = samples.sparse_variable_relative_weights(i_samples,:);
elseif isprop(samples,'sparse_variable_relative_weights') && isstruct(samples.sparse_variable_relative_weights)
    samples.sparse_variable_relative_weights.(field) = samples.sparse_variable_relative_weights.(field)(i_samples,:);
end

E = NaN(1,size(x,2));

for d=1:n_dims
    if issparse(x)
        x_local = nonzeros(x(:,d));
    else
        x_local = x(:,d);
    end
    if isfield(samples.var,'num_instances')
        % Support of outdated form
        [i,~,counts] = find(samples.var.num_instances(:,d));
        w_local = counts.*w_particles(i);
        w_local = w_local/sum(w_local);
    elseif isfield(samples.var,'relative_weights') || isstruct(samples) || isempty(samples.sparse_variable_relative_weights)
        % Includes support of outdated form
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