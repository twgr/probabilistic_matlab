classdef mv_gaussian_class < base_primitive
    properties
        mu      % Always an array of size NxM
        sigma   % Can either be one array for all possible mu, or a MxMxN
        % array.  Vector operation possible when only single sigma,
        % otherwise need to loop over the sigma's to calculate the
        % cholesky decomposition which will be slow.
        Z       % Optional input, allows truncation of the distribution in
        % the cumulative normal space. Of size Mx2xN
    end
    
    properties (Hidden=true, SetAccess=private)
        chol_sigma
        chol_sigma_inv
        det_sigma
        log_norm_constant
    end
    
    methods
        
        function obj = mv_gaussian_class(mu,sigma,Z)
            obj.mu      = mu;
            if ndims(sigma)==3 && size(sigma,1) ~= size(sigma,2)
                % The sigma parameter is MxMxN for consistency when there
                % it is fixed rather than varying by sample.  However, it
                % is natural to provide it as NxMxM.  This autocorrects
                % when it is provided like this
                sigma = permute(sigma,[3,2,1]);
            end
            obj.sigma   = sigma;
            if exist('Z','var')
                obj.Z = Z;
            end
        end
        
        function obj = set.sigma(obj,sigma)
            if ndims(sigma)==3 && size(sigma,1) ~= size(sigma,2)
                % The sigma parameter is MxMxN for consistency when there
                % it is fixed rather than varying by sample.  However, it
                % is natural to provide it as NxMxM.  This autocorrects
                % when it is provided like this
                sigma = permute(sigma,[3,2,1]);
            end
            obj.sigma = sigma;
            if ismatrix(sigma)
                obj.chol_sigma = chol(sigma); %#ok<*MCSUP>
                obj.chol_sigma_inv = obj.chol_sigma\eye(size(sigma,1));
                obj.det_sigma = (prod(diag(obj.chol_sigma)))^2;
            else
                obj.chol_sigma = NaN(size(sigma));
                obj.chol_sigma_inv = NaN(size(sigma));
                obj.det_sigma = NaN(size(sigma,1),1);
                for n=1:size(sigma,3)
                    obj.chol_sigma(:,:,n) = chol(squeeze(sigma(:,:,n)));
                    obj.chol_sigma_inv(:,:,n) = squeeze(obj.chol_sigma(:,:,n))\eye(size(sigma,1));
                    obj.det_sigma(n) = (prod(diag(obj.chol_sigma(:,:,n))))^2;
                end
            end
            if isempty(obj.Z)
                log_Z_trun = 0;
            else
                log_Z_trun = sum(log(squeeze(abs(obj.Z(:,2,:)-obj.Z(:,1,:)))),1)';
            end
            obj.log_norm_constant = log_Z_trun+0.5*size(obj.mu,2)*log(2*pi)+0.5*log(obj.det_sigma);
        end
        
        function obj = set.Z(obj,Z)
            if ndims(Z)==3 && size(Z,1)~=size(obj.sigma,1)
                % Z provided the wrong way around so flip.
                Z = permute(Z,[3,2,1]);
            elseif ndims(Z)==3 && size(Z,1)==size(obj.sigma,1) && size(Z,3)==size(obj.sigma,1)
                % Not clear if correct or not.
                warning('Z is Nx2xN, not clear which orientation provided.  Assumes 1st dim is sample id by default');
                Z = permute(Z,[3,2,1]);
            end                
            obj.Z = Z;
            if isempty(obj.Z)
                log_Z_trun = 0;
            else
                log_Z_trun = sum(log(squeeze(abs(obj.Z(:,2,:)-obj.Z(:,1,:)))),1)';
            end
            obj.log_norm_constant = log_Z_trun+0.5*size(obj.mu,2)*log(2*pi)+0.5*log(obj.det_sigma);
        end
        
        function vals = draw(obj,n_draws)
            assert(any(size(obj.mu,1)==[1,n_draws]) && any(size(obj.sigma,3)==[1,n_draws]),...
                'Obj must either have single value for parameters or the same number as wish to be sampled');
            D = size(obj.mu,2);
            
            if isempty(obj.Z)
                phis = randn(n_draws,size(obj.mu,2));
            else
                Zs = rand(n_draws,D).*squeeze(obj.Z(:,2,:)-obj.Z(:,1,:))'+squeeze(obj.Z(:,1,:))';
                phis = norminv(Zs,0,1);
            end
            if ismatrix(obj.chol_sigma)
                vals = bsxfun(@plus,obj.mu,phis*obj.chol_sigma);
            else
                vals = bsxfun(@plus,obj.mu,...
                    squeeze(sum(bsxfun(@times,phis,...
                    permute(obj.chol_sigma,[3,1,2])),2)));
            end
        end
        
        function log_p = log_pdf(obj,vals)
            d = bsxfun(@minus,vals,obj.mu);
            if ismatrix(obj.chol_sigma_inv)
                s = d*obj.chol_sigma_inv;
            else
                s = squeeze(sum(bsxfun(@times,d,permute(obj.chol_sigma_inv,[3,1,2])),2));
            end
            log_p = -obj.log_norm_constant-0.5*(sum(s.^2,2));
        end
        
        
        function p = pdf(obj,vals)
            p = exp(obj.observe(vals));
        end
        
        function c = cdf(obj,vals)
            % By convention this is the cdf of the marginal of the first
            % variable then the cdf of the second conditioned on the first
            % and so on.  As such it returns columns equal to the
            % dimensionality.
            d = bsxfun(@minus,vals,obj.mu);
            if ismatrix(obj.chol_sigma_inv)
                s = d*obj.chol_sigma_inv;
            else
                s = squeeze(sum(bsxfun(@times,d,permute(obj.chol_sigma_inv,[3,1,2])),2));
            end
            c = normcdf(s);
        end
        
        function log_c = log_cdf(obj,vals)
            log_c = log(obj.cdf(vals));
        end
    end
end