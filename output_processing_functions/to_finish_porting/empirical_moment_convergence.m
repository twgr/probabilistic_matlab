function varargout = empirical_moment_convergence(samples,orders,bIgnoreNaN,varargin)
%varargout = empirical_moment_convergence(samples,order,bIgnoreNaN,varargin)
%
% Calculations the moments of order orders the variables named in vararign 
% within the sample set samples.  Outputs are grouped first by the same
% varargin such that
%   empirical_moment(samples,[1,3],true,'x',y')
% returns [mean_x, skewness_x, mean_y, skewness_y]
%
% order = 1 : mean
% order = 2 : standard deviation
% order = 3 : skewness
% order = 4 : excess kurtosis
% order = n (>4) : nth order moment given by
%                       E[(X-muX).^n]./(std_dev.^n) with appropriate
%                       weighting of the samples incorporated.


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
        varargout{n} = empirical_mean(samples,[],bIgnoreNaN,varargin{n});
    end
    return
end

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
iteration_reorder = iteration_reorder(:);

if isempty(samples.sparse_variable_relative_weights)
    % This is the current case without compression
    n_out = 1;
    
    for n=1:numel(varargin)        
        X = samples.var.(varargin{n})(iteration_reorder,:);        
        
        if bIgnoreNaN
            bNaN = isnan(X);
            X(bNaN) = 0;
            w = samples.relative_particle_weights(iteration_reorder,:);
            w =  bsxfun(@rdivide,w,sum(bsxfun(@times,~bNaN,w),1));
            w(bNaN) = 0;
        else
            w = samples.relative_particle_weights(iteration_reorder,:);
            w = w/sum(w);
        end       
        scales = cumsum(w,1);
        scales = scales(n_iter:n_iter:end,:);
        cumX = cumsum(bsxfun(@times,X,w),1);
        meanX = cumX(n_iter:n_iter:end,:)./scales;
        
        E_X2 = cumsum(bsxfun(@times,X.^2,w),1);
        E_X2 = E_X2(n_iter:n_iter:end,:)./scales;
        
        varX = E_X2-meanX.^2;
        varX = max(0,varX);
        
        E_X = cell(max(orders),1);
        for no = 3:max(orders)
            E_X{no} = cumsum(bsxfun(@times,X.^no,w),1);
            E_X{no} = E_X{no}(n_iter:n_iter:end,:)./scales;
        end
        
        for no = 1:numel(orders)
            if orders(no)==1
                varargout{n_out} = meanX;
            elseif orders(no)==2
                varargout{n_out} = sqrt(varX);
            else
                varargout{n_out} = (-1).^(no-1).*(no-1)*meanX.^no+(-1).^(no).*nchoosek(no,no-2).*E_X2.*meanX.^(no-2);
                for mu_power = 0:(no-3)
                    varargout{n_out} = varargout{n_out}+nchoosek(no,mu_power).*((-meanX).^mu_power).*E_X{no-mu_power};
                end
                varargout{n_out} = varargout{n_out}./(varX.^(orders(no)/2));
                if orders(no)==4
                    varargout{n_out} = varargout{n_out} - 3;
                end
            end
            varargout{n_out} = full(varargout{n_out});
            n_out = n_out+1;
        end
    end
elseif isnumeric(samples.sparse_variable_relative_weights)
    % Current case with compression but common sparsity structure
    w = samples.sparse_variable_relative_weights(iteration_reorder,:);
    sw = size(w);
    iNonZeros = find(w);
    w = full(w(iNonZeros));    
    [iM,jM] = ind2sub(sw,iNonZeros);
    
    n_out = 1;
    for n=1:numel(varargin)        
        X = samples.var.(varargin{n})(iteration_reorder,:);
        [~,sX2] = size(X);
        X = full(X(iNonZeros));
        w_local = w;
        iM_local = iM;
        jM_local = jM;
        if bIgnoreNaN
            bNotNaN = ~isnan(X);
            w_local = w(bNotNaN);
            X = X(bNotNaN);
            iM_local = iM(bNotNaN);
            jM_local = jM_local(bNotNaN);
        end
        i_iter = ceil(iM_local/(n_particles*n_nodes));        
        scales = cumsum(accumarray([i_iter, jM_local],w_local,[n_iter,sX2]),1);
        meanX = cumsum(accumarray([i_iter, jM_local],w_local.*X,[n_iter,sX2]))./scales;
        
        E_X2 = cumsum(accumarray([i_iter, jM_local],w_local.*X.^2,[n_iter,sX2]))./scales;
        varX = E_X2-meanX.^2;
        varX = max(0,varX);
        
        E_X = cell(max(orders),1);
        for no = 3:max(orders)
            E_X{no} = cumsum(accumarray([i_iter, jM_local],w_local.*X.^no,[n_iter,sX2]))./scales;
        end
        
        for no = 1:numel(orders)
            if orders(no)==1
                varargout{n_out} = meanX;
            elseif orders(no)==2
                varargout{n_out} = sqrt(varX);
            else
                varargout{n_out} = (-1).^(no-1).*(no-1)*meanX.^no+(-1).^(no).*nchoosek(no,no-2).*E_X2.*meanX.^(no-2);
                for mu_power = 0:(no-3)
                    varargout{n_out} = varargout{n_out}+nchoosek(no,mu_power).*((-meanX).^mu_power).*E_X{no-mu_power};
                end
                varargout{n_out} = varargout{n_out}./(varX.^(orders(no)/2));
                if orders(no)==4
                    varargout{n_out} = varargout{n_out} - 3;
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
        
        X = samples.var.(varargin{n})(iteration_reorder,:);
        [sX1,sX2] = size(X);
        
        w_local = samples.sparse_variable_relative_weights.(varargin{n})(iteration_reorder,:);
        iNonZeros = find(w_local);
        w_local = w_local(iNonZeros);
        X = X(iNonZeros);
        [iM_local,jM_local] = ind2sub([sX1,sX2],iNonZeros);        
        if bIgnoreNaN
            bNotNaN = ~isnan(X);
            w_local = w_local(bNotNaN);
            X = X(bNotNaN);
            iM_local = iM_local(bNotNaN);
            jM_local = jM_local(bNotNaN);
        end
        
        i_iter = ceil(iM_local/(n_particles*n_nodes));        
        scales = cumsum(accumarray([i_iter, jM_local],w_local,[n_iter,sX2]),1);
        meanX = cumsum(accumarray([i_iter, jM_local],w_local.*X,[n_iter,sX2]))./scales;
        
        E_X2 = cumsum(accumarray([i_iter, jM_local],w_local.*X.^2,[n_iter,sX2]))./scales;
        varX = E_X2-meanX.^2;
        varX = max(0,varX);
        
        E_X = cell(max(orders),1);
        for no = 3:max(orders)
            E_X{no} = cumsum(accumarray([i_iter, jM_local],w_local.*X.^no,[n_iter,sX2]))./scales;
        end
        
        for no = 1:numel(orders)
            if orders(no)==1
                varargout{n_out} = meanX;
            elseif orders(no)==2
                varargout{n_out} = sqrt(varX);
            else
                varargout{n_out} = (-1).^(no-1).*(no-1)*meanX.^no+(-1).^(no).*nchoosek(no,no-2).*E_X2.*meanX.^(no-2);
                for mu_power = 0:(no-3)
                    varargout{n_out} = varargout{n_out}+nchoosek(no,mu_power).*((-meanX).^mu_power).*E_X{no-mu_power};
                end
                varargout{n_out} = varargout{n_out}./(varX.^(orders(no)/2));
                if orders(no)==4
                    varargout{n_out} = varargout{n_out} - 3;
                end
            end
            varargout{n_out} = full(varargout{n_out});
            n_out = n_out+1;
        end
    end
end
    