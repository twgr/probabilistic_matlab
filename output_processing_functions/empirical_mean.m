function varargout = empirical_mean(samples,i_samples,bIgnoreNaN,varargin)
%empirical_mean
%
% varargout = empirical_mean(samples,i_samples,bIgnoreNaN,varargin)
% 
% Calculates the empirical average of variables in a stack_object.
% Averages are taken seperately across each dimension of the specified
% variables.
%
% Inputs:
%   samples = stack_object to take empirical mean from
%   i_samples = Subset of particles to take in the average.  If left blank
%               then all are used.
%   bIgnoreNaN = If true NaNs will be ignored in the averaging, if false
%                then the mean will be NaN if there is any NaN in the
%                samples.  If left empty, defaults to false.
%   varargin = Strings giving names of variables to take empirical mean
%              for. One output is provided for each provided string. There
%              are computational savings to making a single call to the
%              function with many variables specified then making numerous
%              calls.
% Outputs:
%   varargout = Empirical means.
%
%
% Example
%   [m_x,m_y] = emprical_mean(samples,[],false,'x','y');
%   
% Tom Rainforth 05/07/16

if isempty(i_samples)
    i_samples = (1:size(samples.var.(varargin{end}),1))';
end

if isempty(bIgnoreNaN)
    bIgnoreNaN = false;
end

if ~all(cellfun(@(x) isnumeric(samples.var.(x)), varargin))
    error('empirical_mean can only be called for numeric variables');
end

varargout = cell(1,numel(varargin));

if isempty(samples.sparse_variable_relative_weights)
    % No compression
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
    % Compression but common sparsity structure
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
    
