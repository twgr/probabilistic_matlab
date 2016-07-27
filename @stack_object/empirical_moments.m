function varargout = empirical_moments(samples,orders,i_samples,bIgnoreNaN,varargin)
%varargout = empirical_moments(samples,order,i_samples,bIgnoreNaN,varargin)
%
% Calculations the moments of order orders the variables named in vararign 
% within the sample set samples.  
%
% Inputs: samples, orders, i_samples, bIgnoreNaN, varargin.
%   All as per empirical_mean other than the order input:
% order = 1 : mean.  As per empirical_mean
% order = 2 : standard deviation
% order = 3 : skewness
% order = 4 : excess kurtosis
% order = n (>4) : nth order moment given by
%                       E[(X-muX).^n]./(std_dev.^n) with appropriate
%                       weighting of the samples incorporated.
%
% Outputs: as per empriical_mean.  Results are grouped first by the same
% varargin such that
%   empirical_moment(samples,[1,3],[],true,'x',y')
% returns [mean_x, skewness_x, mean_y, skewness_y]
%
% Tom Rainforth 05/07/16

if ~exist('i_samples','var') || isempty(i_samples) || ischar(i_samples) || numel(i_samples)==1    
    if ischar(i_samples)
        if exist('bIgnoreNaN','var')
            varargin = [{i_samples},{bIgnoreNaN}, varargin];
        else
            varargin = {i_samples};
        end
        bIgnoreNaN = false;        
    elseif numel(i_samples)==1
        varargin = [{bIgnoreNaN}, varargin];
        bIgnoreNaN = i_samples;
    end
    if ~isempty(varargin)
        i_samples = (1:size(samples.var.(varargin{end}),1))';
    else
        i_samples = (1:size(samples.var.(bIgnoreNaN),1))';
    end
end

if ~exist('bIgnoreNaN','var') || isempty(bIgnoreNaN)
    bIgnoreNaN = false;
elseif ischar(bIgnoreNaN)
    varargin = [{bIgnoreNaN}, varargin];
    bIgnoreNaN = false;
end

if ~all(cellfun(@(x) isnumeric(samples.var.(x)), varargin))
    error('Only valid for numeric variables');
end

varargout = cell(1,numel(varargin));

if numel(orders)==1 && orders==1
    for n=1:numel(varargin)
        varargout{n} = empirical_mean(samples,i_samples,bIgnoreNaN,varargin{n});
    end
    return
end


if isempty(samples.sparse_variable_relative_weights)
    % No compression
    n_out = 1;
    
    for n=1:numel(varargin)        
        X = samples.var.(varargin{n})(i_samples,:);        
        if bIgnoreNaN
            bNaN = isnan(X);
            X(bNaN) = 0;
            w =  bsxfun(@rdivide,samples.relative_particle_weights(i_samples,:),sum(bsxfun(@times,~bNaN,samples.relative_particle_weights(i_samples,:)),1));
            w(bNaN) = 0;
        else
            w = samples.relative_particle_weights(i_samples,:)/sum(samples.relative_particle_weights(i_samples,:));
        end
        scale = sum(w,1);
        meanX = sum(bsxfun(@times,X,w),1)./scale;
        X = bsxfun(@minus,X,meanX);
        varX = sum(bsxfun(@times,w,X.^2),1)./scale;
        varX = max(0,varX);
        for no = 1:numel(orders)     
            if orders(no)==1
                varargout{n_out} = meanX;
            elseif orders(no)==2
                varargout{n_out} = sqrt(varX);
            else
                varargout{n_out} = sum(bsxfun(@times,w,X.^orders(no)),1)./(scale.*varX.^(orders(no)/2));
                if orders(no)==4
                    varargout{n_out} = varargout{n_out} - 3;
                end
            end
            varargout{n_out} = full(varargout{n_out});
            n_out = n_out+1;
        end
    end
elseif isnumeric(samples.sparse_variable_relative_weights)
    % Compression but common sparsity structure
    w = samples.sparse_variable_relative_weights(i_samples,:);
    sw = size(w);
    iNonZeros = find(w);
    w = full(w(iNonZeros));    
    [~,jM] = ind2sub(sw,iNonZeros);
    
    n_out = 1;
    for n=1:numel(varargin)        
        X = samples.var.(varargin{n})(i_samples,:);
        [~,sX2] = size(X);
        X = X(iNonZeros);
        w_local = w;
        jM_local = jM;
        if bIgnoreNaN
            bNotNaN = ~isnan(X);
            w_local = w(bNotNaN);
            X = X(bNotNaN);
            jM_local = jM_local(bNotNaN);
        end
        scale = accumarray(jM_local,w_local,[sX2,1]);
        meanX = accumarray(jM_local,w_local.*X,[sX2,1])./scale;
        X = X-meanX(jM_local);
        varX = accumarray(jM_local,w_local.*(X.^2),[sX2,1])./scale;
        varX = max(0,varX);
        for no = 1:numel(orders)
            if orders(no)==1
                varargout{n_out} = meanX;
            elseif orders(no)==2
                varargout{n_out} = sqrt(varX);
            else
                varargout{n_out} = accumarray(jM_local,w_local.*(X.^orders(no)),[sX2,1])./(scale.*varX.^(orders(no)/2));
                if orders(no)==4
                    varargout{n_out} = varargout{n_out}-3;
                end
            end
            varargout{n_out} = full(varargout{n_out});
            n_out = n_out+1;
        end
    end
else
    % Different variables have a different sparsity structure
    n_out = 1;
    
    for n=1:numel(varargin)
        
        X = samples.var.(varargin{n})(i_samples,:);
        sX = size(X);
        
        w_local = samples.sparse_variable_relative_weights.(varargin{n})(i_samples,:);
        iNonZeros = find(w_local);
        w_local = w_local(iNonZeros);
        X = X(iNonZeros);
        [~,jM_local] = ind2sub(sX,iNonZeros);        
        if bIgnoreNaN
            bNotNaN = ~isnan(X);
            w_local = w_local(bNotNaN);
            X = X(bNotNaN);
            jM_local = jM_local(bNotNaN);
        end
        
        scale = accumarray(jM_local,w_local,[sX2,1]);
        meanX = accumarray(jM_local,w_local.*X,[sX2,1])./scale;
        X = X-meanX(jM_local);
        varX = accumarray(jM_local,w_local.*(X.^2),[sX2,1])./scale;
        varX = max(0,varX);
        for no = 1:numel(orders)
            if orders(no)==1
                varargout{n_out} = meanX;
            elseif orders(no)==2
                varargout{n_out} = sqrt(varX);
            else
                varargout{n_out} = accumarray(jM_local,w_local.*(X.^orders(no)),[sX2,1])./(scale.*varX.^(orders(no)/2));
                if orders(no)==4
                    varargout{n_out} = varargout{n_out}-3;
                end
            end
            varargout{n_out} = full(varargout{n_out});
            n_out = n_out+1;
        end
    end
end
    
