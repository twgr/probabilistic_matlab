function [densities,edges] = results_binning(samples,field,d_start,d_end,n_bins_or_edges,b_discrete)

densities = cell(1+d_end-d_start,1);

if isscalar(n_bins_or_edges)
    n_bins = n_bins_or_edges;
    edges = cell(1+d_end-d_start,1);
else
    if b_discrete
        n_bins = numel(n_bins_or_edges);
    else
        n_bins = numel(n_bins_or_edges)-1;
    end
    edges = repmat({n_bins_or_edges},1+d_end-d_start,1);
end

for d=d_start:d_end
    if issparse(samples.var.(field))
        X = nonzeros(samples.var.(field)(:,d));
    else
        X = samples.var.(field)(:,d);
    end
    
    if isscalar(n_bins_or_edges)
        maxX = max(X);
        minX = min(X);
        if ~b_discrete
            maxX = minX+(1/n_bins_or_edges)*(maxX-minX);
            edges{d}=(linspace(minX,maxX,n_bins_or_edges+1));
        else
            edges{d}=(linspace(minX,maxX,n_bins_or_edges));
        end        
    else
        minX = min(edges{d}(1),edges{d}(end));
        maxX = max(edges{d}(1),edges{d}(end));
    end
        
    if isfield(samples.var,'num_instances')
        % Support of outdated form
        [i,~,counts] = find(samples.var.num_instances(:,d));
        w = counts.*samples.relative_weights(i);
    elseif isfield(samples.var,'relative_weights') || isstruct(samples)
        % Support of outdated form
        w = samples.relative_weights;
    elseif isempty(samples.sparse_variable_relative_weights)
        if isempty(samples.relative_particle_weights)
            w = ones(size(X,1),1);
        else
            w = samples.relative_particle_weights;
        end
    elseif isnumeric(samples.sparse_variable_relative_weights)
        w = nonzeros(samples.sparse_variable_relative_weights(:,d));
    else
        w = nonzeros(samples.sparse_variable_relative_weights.(field)(:,d));
    end
    
    w = w/sum(w);
    
    [~, bin_assignment] = histc(full(X),edges{d});
    counts = accumarray(bin_assignment,w);
    densities{d} = zeros(n_bins,1);
    den_incomp = counts*numel(counts)/(maxX-minX);
    densities{d}(1:(min(n_bins,numel(den_incomp)))) = den_incomp(1:(min(n_bins,numel(den_incomp))));
end

densities = [densities{:}]';

if exist('b_discrete','var') && b_discrete
    densities = bsxfun(@rdivide,densities,sum(densities,2));
end