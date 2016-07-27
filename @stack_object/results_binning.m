function [densities,edges,b_distinct_values] = results_binning(samples,field,n_bins,b_discrete,dims,i_samples,edges)
%results_binning Seperates samples into bins
%
% [densities,edges,b_distinct_values] = ...
%   results_binning(samples,field,d_start,d_end,n_bins_or_edges,b_discrete)
%
% Inputs:
%   samples = stack_object
%   field = Field to bin
%   n_bins = Number of bins to use.  Can be left empty if edges is set
%   b_discrete = Whether to treat the inputs are discrete variables at
%                return probabilities instead of densities (i.e. take a
%                counting measure instead of Lebesgue).  Default = false;
%   dims = Dimensions to calculate the ess for.  By default, all are
%          calculated.
%   i_samples = The indices of the samples to use in the calculate.  By
%               default, all are used.
%   edges = If not set, n_bins evenly spaced bins are used.  Otherwise the
%           bin edges can be manually set.  Can even be a vector of edges
%           (as per the histc command) or a cell array of edges, defining
%           different edges for each of the dimensions.  If b_discrete =
%           true, then the edges are the discrete points to take counts at
%           rather than the edges.
%
% Outputs:
%   densities = Density estimates if b_discrete=false, otherwise
%               probability estimates.
%   edges = As per the input (allows returning if these were internally
%           derived).
%   b_distinct_values = List of booleans indicating whether there was more
%                       than one unique value along that dimension.
%
%  Tom Rainforth 27/07/16

[n_samples,n_dims_total] = size(samples.var.(field));

if ~exist('b_discrete','var') || isempty(b_discrete)
    b_discrete = false;
else
    assert(isscalar(b_discrete),'b_discrete should be a scalar');
end

if ~exist('dims','var') || isempty(dims)
    dims = 1:n_dims_total;
end
ndims = numel(dims);

if ~exist('i_samples','var') || isempty(i_samples)
    i_samples = (1:n_samples)';
end

if ~exist('edges','var') || isempty(edges)
    edges = cell(ndims,1);
elseif isnumeric(edges)
    edges = repmat({edges},ndims,1);
elseif ~isempty(nbins) && nbins~=2
    for n=1:numel(edges)
        if numel(edges{n})==2
            edges{n} = linspace(min(edges{n}),max(edges{n}),nbins)';
        end
    end
end

densities = cell(ndims,1);
b_distinct_values = true(ndims,1);

for nd = 1:numel(dims)
    d=dims(nd);
    
    x_local = samples.get_variable_dim(field,d,i_samples);
    
    if isempty(edges{nd})
        maxX = max(x_local);
        minX = min(x_local);
        if maxX==minX
            % Only 1 unique value
            edges{nd} = x_local(1)*ones(n_bins,1);
            if b_discrete
                densities = [1;zeros(n_bins-1,1)];
            else
                densities = NaN(size(edges{nd}));
            end
            b_distinct_values(nd) = false;
            continue
        elseif ~b_discrete
            maxX = minX+((n_bins+1e-8)/n_bins)*(maxX-minX);
            edges{nd}=(linspace(minX,maxX,n_bins+1))';
        else
            edges{nd}=(linspace(minX,maxX,n_bins))';
        end
        
    else
        minX = min(edges{nd}(1),edges{nd}(end));
        maxX = max(edges{nd}(1),edges{nd}(end));
        n_bins = numel(edges{nd});
        if ~b_discrete
            n_bins = n_bins-1;
        end
    end
    edges{nd} = edges{nd}(:);    
    
    w_local = samples.get_weights(field,d,i_samples);
    
    [~, bin_assignment] = histc(full(x_local),edges{nd});
    counts = accumarray(bin_assignment,w_local);
    densities{nd} = zeros(n_bins,1);
    den_incomp = counts*numel(counts)/(maxX-minX);
    densities{nd}(1:(min(n_bins,numel(den_incomp)))) = den_incomp(1:(min(n_bins,numel(den_incomp))));
end

densities = [densities{:}];
edges = [edges{:}];

if b_discrete
    densities = bsxfun(@rdivide,densities,sum(densities,1));
end