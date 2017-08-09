function [densities,edges,b_distinct_values,counts,variances] = results_binning(samples,field,n_bins,b_discrete,dims,i_samples,edges,bf)
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
%           rather than the edges.  If edges are a matrix, the first
%           dimension should be the different edges of one binning
%           strategy.
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
n_dims = numel(dims);

if ~exist('i_samples','var') || isempty(i_samples)
    i_samples = (1:n_samples)';
end

if ~exist('edges','var') || isempty(edges)
    edges = cell(n_dims,1);
elseif isnumeric(edges)
    edges = repmat({edges},n_dims,1);
elseif ~isempty(n_bins) && n_bins~=2
    for n=1:numel(edges)
        if numel(edges{n})==2
            edges{n} = linspace(min(edges{n}),max(edges{n}),n_bins)';
        end
    end
end

densities = cell(n_dims,1);
b_distinct_values = true(n_dims,1);
bMatInput = false;
bin_assignments = cell(n_dims,1);
variances = cell(n_dims,1);
counts = cell(n_dims,1);

for nd = 1:numel(dims)
    d=dims(nd);
    
    x_local = samples.get_variable(field,d,i_samples);
    
    if isempty(edges{nd})
        maxX = max(x_local);
        minX = min(x_local);
        if maxX==minX
            % Only 1 unique value
            edges{nd} = x_local(1)*ones(n_bins,1);
            if b_discrete
                densities{nd} = [1;zeros(n_bins-1,1)];
            else
                densities{nd} = NaN(size(edges{nd}));
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
        bMatInput = bMatInput || ~isvector(edges{nd});
        
        if size(edges{nd},1)==1
            edges{nd} = edges{nd}(:);
        end
        %
        %         minX = min(edges{nd},[],1);
        %         maxX = max(edges{nd},[],1);
        n_bins = size(edges{nd},1);
        if ~b_discrete
            n_bins = n_bins-1;
        end
    end
    
    w_local = samples.get_weights(field,d,i_samples);
    x_local = full(x_local);
    if exist('bf','var') && bf
        f_local = samples.get_variable('fs',1,i_samples);
        fw_local = samples.get_weights('fs',1,i_samples);
    end
    bin_assignments{nd} = histmat(x_local,edges{nd});
    n_ass = size(bin_assignments{nd},2);
    counts{nd} = zeros(n_bins,n_ass);
    variances{nd} = zeros(n_bins,n_ass);
    for n=1:size(counts{nd},2)
        i_in = bin_assignments{nd}(:,n)~=0;
        counts{nd}(:,n) = accumarray(bin_assignments{nd}(i_in,n),w_local(i_in),[n_bins,1]);
        if exist('bf','var') && bf
            f1 = accumarray(bin_assignments{nd}(i_in,n),fw_local(i_in).*f_local(i_in),[n_bins,1]);    
            f2 = accumarray(bin_assignments{nd}(i_in,n),fw_local(i_in).*f_local(i_in).^2,[n_bins,1]);        
            V1 = accumarray(bin_assignments{nd}(i_in,n),fw_local(i_in),[n_bins,1]);    
            V2 = accumarray(bin_assignments{nd}(i_in,n),fw_local(i_in).^2,[n_bins,1]);  
            variances{nd}(:,n) = f2./(V1-V2./V1)-f1./(V1.^2-V1);
            variances{nd}(isnan(variances{nd}(:,n)),n)=0;
        else
            variances{nd}(:,n) = var(x_local(i_in));
            variances{nd}(isnan(variances{nd}(:,n)),n)=0;
        end
    end
    if b_discrete
        densities{nd} = bsxfun(@rdivide,counts{nd},sum(counts{nd},1));
    else
        densities{nd} = bsxfun(@rdivide,counts{nd},diff(edges{nd},[],1));
    end
    %    densities{nd} = bsxfun(@rdivide,counts*n_bins,(maxX-minX));
    %     densities{nd} = zeros(n_bins,n_ass);
    %     den_incomp = bsxfun(@rdivide,counts*n_bins,(maxX-minX));
    %densities{nd}(1:(min(n_bins,numel(den_incomp)))) = den_incomp(1:(min(n_bins,numel(den_incomp))));
end

if ~bMatInput
    densities = [densities{:}];
    edges = [edges{:}];
    counts = [counts{:}];
    variances = [variances{:}];
end



end

function bin_assignment = histmat(x,edges)

if ~isvector(edges)
    b_greater = bsxfun(@ge,x,reshape(edges,1,size(edges,1),[]));
    bin_assignment = sum(b_greater,2);
    if size(bin_assignment,1)==1
        bin_assignment = squeeze(bin_assignment)';
    else
        bin_assignment = squeeze(bin_assignment);
    end
else
    [~, bin_assignment] = histc(x,edges);
end

end