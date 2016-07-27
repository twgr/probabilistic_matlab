function [densities,edges] = results_binning(samples,field,n_bins,b_discrete,dims,i_samples,edges)
%results_binning Seperates samples into bins
%
% [densities,edges] = results_binning(samples,field,d_start,...
%                                         d_end,n_bins_or_edges,b_discrete)
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
%
%  Tom Rainforth 27/07/16

x = samples.var.(field);
n_samples_total = size(x,1);

if ~exist('b_discrete','var') || isempty(b_discrete)
    b_discrete = false;
end

if ~exist('dims','var') || isempty(dims)
    dims = 1:size(x,2);
end

if ~exist('i_samples','var') || isempty(i_samples)
    i_samples = (1:n_samples_total)';
end

ndims = numel(dims);

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

if ~exist('edges','var') || isempty(edges)
    edges = cell(ndims,1);
elseif isnumeric(edges)
    edges = repmat({edges},ndims,1);
end

densities = cell(ndims,1);

for d=dims(:)'
    if issparse(x)
        x_local = nonzeros(x(:,d));
    else
        x_local = x(:,d);
    end
    
    if isempty(edges{d})
        maxX = max(x_local);
        minX = min(x_local);
        if ~b_discrete            
            maxX = minX+((n_bins+1e-8)/n_bins)*(maxX-minX);
            edges{d}=(linspace(minX,maxX,n_bins+1));
        else
            edges{d}=(linspace(minX,maxX,n_bins));
        end        
    else
        minX = min(edges{d}(1),edges{d}(end));
        maxX = max(edges{d}(1),edges{d}(end));
        n_bins = numel(edges{d});
        if ~b_discrete
            n_bins = n_bins-1;
        end
    end
    
    if isempty(samples.sparse_variable_relative_weights)
        w_local = w_particles;
    elseif isnumeric(samples.sparse_variable_relative_weights)
        w_local = nonzeros(samples.sparse_variable_relative_weights(:,d));
    else
        w_local = nonzeros(samples.sparse_variable_relative_weights.(field)(:,d));
    end
    
    w_local = w_local/sum(w_local);
    
    [~, bin_assignment] = histc(full(x_local),edges{d});
    counts = accumarray(bin_assignment,w_local);
    densities{d} = zeros(n_bins,1);
    den_incomp = counts*numel(counts)/(maxX-minX);
    densities{d}(1:(min(n_bins,numel(den_incomp)))) = den_incomp(1:(min(n_bins,numel(den_incomp))));
end

densities = [densities{:}]';

if b_discrete
    densities = bsxfun(@rdivide,densities,sum(densities,2));
end