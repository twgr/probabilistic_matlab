function varargout = empirical_moments_convergence(samples,orders,field,dims,bIgnoreNaN)
%varargout = empirical_moment_convergence(samples,orders,field,dims,bIgnoreNaN)
%
% As empirical_moments except calculates the values for each MCMC iteration
% and does not allow the i_samples input (as this would severly complicate
% things and is not really necessary).
%
% Tom Rainforth 05/07/16


[n_samples,n_dims_total] = size(samples.var.(field));

if ~exist('dims','var') || isempty(dims)
    dims = 1:n_dims_total;
end

i_samples = (1:n_samples)';

if ~exist('bIgnoreNaN','var') || isempty(bIgnoreNaN)
    bIgnoreNaN = false;
end

if ~isnumeric(samples.var.(field))
    error('Only valid for numeric variables');
end

n_iter = samples.options.n_iter;
n_p_per_iter = n_samples/n_iter;

varargout = cell(1,numel(orders));

X = samples.get_variable(field,dims,i_samples);
[w,~,iNonZeros] = samples.get_weights(field,dims,i_samples);
[iM,jM] = ind2sub([n_samples,n_dims_total],iNonZeros);
if bIgnoreNaN
    bNotNaN = ~isnan(X);
    w = w(bNotNaN);
    X = X(bNotNaN);
    iM = iM(bNotNaN);
    jM = jM(bNotNaN);
end
i_iter = ceil(iM/(n_p_per_iter));  
scales = cumsum(accumarray([i_iter, jM],w,[n_iter,n_dims_total]),1);
meanX = cumsum(accumarray([i_iter, jM],w.*X,[n_iter,n_dims_total]))./scales;
if numel(orders)==1 && orders==1
    varargout{1} = meanX;
    return
end

% In this convergence case, need to use a more complicated formula then the
% standard (or will be crazy inefficient as have do a for loop over
% iterations).  The below is based on expanding out the binominal series
% inside of E[(X-mu_X)^n] to give a number of terms of E[X], E[X^2],
% E[X^3], etc that can each be calculated in a vectorized fashion.
E_X2 = cumsum(accumarray([i_iter, jM],w.*X.^2,[n_iter,n_dims_total]))./scales;
varX = E_X2-meanX.^2;
varX = max(0,varX);

E_X = cell(max(orders),1);
for no = 3:max(orders)
    E_X{no} = cumsum(accumarray([i_iter, jM],w.*X.^no,[n_iter,n_dims_total]))./scales;
end

for no = 1:numel(orders)
    if orders(no)==1
        varargout{no} = meanX;
    elseif orders(no)==2
        varargout{no} = sqrt(varX);
    else
        varargout{no} = (-1).^(no-1).*(no-1)*meanX.^no+(-1).^(no).*nchoosek(no,no-2).*E_X2.*meanX.^(no-2);
        for mu_power = 0:(no-3)
            varargout{no} = varargout{no}+nchoosek(no,mu_power).*((-meanX).^mu_power).*E_X{no-mu_power};
        end
        varargout{no} = varargout{no}./(varX.^(orders(no)/2));
        if orders(no)==4
            varargout{no} = varargout{no} - 3;
        end
    end
    varargout{no} = full(varargout{no});
end


end
