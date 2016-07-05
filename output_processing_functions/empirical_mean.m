function varargout = empirical_mean(samples,i_samples,bIgnoreNaN,varargin)

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
    i_samples = (1:size(samples.var.(varargin{end}),1))';
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

if isempty(samples.sparse_variable_relative_weights)
    % This is the current case without compression
    for n=1:numel(varargin)
        
        X = samples.var.(varargin{n})(i_samples,:);
        
        if bIgnoreNaN
            bNaN = isnan(X);
            X(bNaN) = 0;
            varargout{n} = sum(bsxfun(@times,samples.relative_particle_weights(i_samples,:),X),1)./...
                sum(bsxfun(@times,~bNaN,samples.relative_particle_weights(i_samples,:)),1);
        else
            varargout{n} = 1/(sum(samples.relative_particle_weights(i_samples,:)))*sum(...
                bsxfun(@times,samples.relative_particle_weights(i_samples,:),X),1);
        end
        varargout{n} = full(varargout{n});
    end
elseif isnumeric(samples.sparse_variable_relative_weights)
    % Current case with compression but common sparsity structure
    w = samples.sparse_variable_relative_weights(i_samples,:);
    sw = size(w);
    iNonZeros = find(w);
    w = full(w(iNonZeros));    
    [~,jM] = ind2sub(sw,iNonZeros);
    for n=1:numel(varargin)        
        X = samples.var.(varargin{n})(i_samples,:);
        XNonZeros = X(iNonZeros);
        if bIgnoreNaN
            bNotNaN = ~isnan(XNonZeros);
            varargout{n} = accumarray(jM,w(bNotNaN).*XNonZeros(bNotNaN),[size(X,2),1])./accumarray(jM,w(bNotNaN),[size(X,2),1]);
        else
            varargout{n} = accumarray(jM,w.*XNonZeros,[size(X,2),1])./accumarray(jM,w,[size(X,2),1]);
        end
        varargout{n} = full(varargout{n});
    end
else
    % Different variables have a different sparsity structure
    for n=1:numel(varargin)
        
        X = samples.var.(varargin{n})(i_samples,:);
        sX = size(X);
        
        w = samples.sparse_variable_relative_weights.(varargin{n})(i_samples,:);
        iNonZeros = find(w);
        w = w(iNonZeros);
        XNonZeros = X(iNonZeros);
        [~,jM] = ind2sub(sX,iNonZeros);
        if bIgnoreNaN
            bNotNaN = ~isnan(XNonZeros);
            varargout{n} = accumarray(jM,w(bNotNaN).*XNonZeros(bNotNaN),[size(X,2),1])./accumarray(jM,w(bNotNaN),[size(X,2),1]);
        else
            varargout{n} = accumarray(jM,w.*XNonZeros,[size(X,2),1])./accumarray(jM,w,[size(X,2),1]);
        end
        varargout{n} = full(varargout{n});
    end
end
    
