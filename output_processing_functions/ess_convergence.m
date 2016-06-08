function E = ess_convergence(samples,field)

x = samples.var.(field);

n_dims = size(x,2);

n_iter = samples.options.n_iter;
if isfield(samples.options,'n_nodes')
    if isfield(samples.options,'b_keep_samples_from_all_nodes') && ~samples.options.b_keep_samples_from_all_nodes
        n_nodes = 1;
    else
        n_nodes = samples.options.n_nodes;
    end
else
    n_nodes = 1;
end

if isfield(samples.options,'b_keep_all_samples_per_node') && ~samples.options.b_keep_all_samples_per_node
    n_particles = 1;
else
    n_particles = samples.options.n_particles;
end
   
iteration_reorder = bsxfun(@plus,(1:n_particles)',(n_particles*n_iter)*(0:1:(n_nodes-1)));
iteration_reorder = bsxfun(@plus,iteration_reorder(:),n_particles*(0:1:(n_iter-1)));

if isfield(samples,'relative_weights') && ~isempty(samples.relative_weights')
    % Support of outdated form
    w_particles = samples.relative_weights(iteration_reorder,:);
    w_particles = w_particles/sum(w_particles);
elseif isprop(samples,'relative_particle_weights')  && ~isempty(samples.relative_particle_weights)
    w_particles = samples.relative_particle_weights(iteration_reorder,:);
else
    w_particles = ones(numel(iteration_reorder),1)/numel(iteration_reorder);
end

x = x(iteration_reorder,:);
if isfield(samples.var,'num_instances')
    samples.var.num_instances = samples.var.num_instances(iteration_reorder,:);
elseif isprop(samples,'sparse_variable_relative_weights') && isnumeric(samples.sparse_variable_relative_weights) && ~isempty(samples.sparse_variable_relative_weights)
    samples.sparse_variable_relative_weights = samples.sparse_variable_relative_weights(iteration_reorder,:);
elseif isprop(samples,'sparse_variable_relative_weights') && isstruct(samples.sparse_variable_relative_weights)
    samples.sparse_variable_relative_weights.(field) = samples.sparse_variable_relative_weights.(field)(iteration_reorder,:);
end

E = NaN(n_iter,size(x,2));


for d=1:n_dims
    if issparse(x)
        iNonZeros = find(x(:,d));
        x_local = nonzeros(x(:,d));
    else
        iNonZeros = 1:size(x_local,2);
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
    
    i_iter = ceil(iNonZeros/(n_particles*n_nodes));  
    w_local = w_local/sum(w_local);
    [~,~,i_val] = fast_unique(x_local);
    V = accumarray([i_iter,i_val],w_local,[n_iter,max(i_val)],[],[],1);
    V = cumsum(V,1);
    E(:,d) = sum(V,2).^2./sum(V.^2,2);
end

end