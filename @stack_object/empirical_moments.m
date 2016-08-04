function varargout = empirical_moments(samples,orders,field,dims,i_samples,bIgnoreNaN)
%varargout = empirical_moments(samples,orders,field,dims,i_samples,bIgnoreNaN)
%
% Calculates empirical moments from the samples
%
% Inputs: samples, orders, field, dims, i_samples, bIgnoreNaN
%   - Order is a vector of the moment orders to calculate where:
%       order = 1 : mean
%       order = 2 : standard deviation
%       order = 3 : skewness
%       order = 4 : excess kurtosis
%       order = n (>4) : nth order moment given by
%                       E[(X-muX).^n]./(std_dev.^n) with appropriate
%                       weighting of the samples incorporated.
%   - For field,dims,i_samples see ess
%   - bIgnoreNaN - If true NaNs will be ignored in the averaging, if false
%     then the mean will be NaN if there is any NaN in the samples.  If
%     left empty, defaults to false.
%
% Outputs: varargout where each output is the empirical moment
% corresponding to a value of order
%
% Tom Rainforth 05/07/16

[n_samples,n_dims_total] = size(samples.var.(field));

if ~exist('dims','var') || isempty(dims)
    dims = 1:n_dims_total;
end

if ~exist('i_samples','var') || isempty(i_samples)
    i_samples = (1:n_samples)';
end

if ~exist('bIgnoreNaN','var') || isempty(bIgnoreNaN)
    bIgnoreNaN = false;
end

if ~isnumeric(samples.var.(field))
    error('Only valid for numeric variables');
end

varargout = cell(1,numel(orders));

X = samples.get_variable_dim(field,dims,i_samples);
[w,~,iNonZeros] = samples.get_weights(field,dims,i_samples);
[~,jM_local] = ind2sub([n_samples,n_dims_total],iNonZeros);
if bIgnoreNaN
    bNotNaN = ~isnan(X);
    w = w(bNotNaN);
    X = X(bNotNaN);
    jM_local = jM_local(bNotNaN);
end
scale = accumarray(jM_local,w,[n_dims_total,1]);
meanX = accumarray(jM_local,w.*X,[n_dims_total,1])./scale;
if numel(orders)==1 && orders==1
    varargout{1} = meanX;
    return
end
X = X-meanX(jM_local);
varX = accumarray(jM_local,w.*(X.^2),[n_dims_total,1])./scale;
varX = max(0,varX);

for no = 1:numel(orders)
    if orders(no)==1
        varargout{no} = meanX;
    elseif orders(no)==2
        varargout{no} = sqrt(varX);
    else
        varargout{no} = accumarray(jM_local,w.*(X.^orders(no)),[n_dims_total,1])./(scale.*varX.^(orders(no)/2));
        if orders(no)==4
            varargout{no} = varargout{no}-3;
        end
    end
    varargout{no} = full(varargout{no});
end

end

