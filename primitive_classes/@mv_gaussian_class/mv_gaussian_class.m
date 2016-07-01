classdef mv_gaussian_class < base_primitive
    properties
        mu      % Always an array of size NxM
        sigma   % Can either be one array for all possible mu, or a MxMxN
                % array.  Vector operation possible when only single sigma,
                % otherwise need to loop over the sigma's to calculate the
                % cholesky decomposition which will be slow.
    end
    
    properties (Hidden=true, SetAccess=private)
        chol_sigma
        chol_sigma_inv
        det_sigma
        log_norm_constant
    end
    
    methods
        
        function obj = mv_gaussian_class(mu,sigma)
            obj.mu      = mu;
            obj.sigma   = sigma;
        end
        
        function obj = set.sigma(obj,sigma)
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
           obj.log_norm_constant = 0.5*size(obj.mu,2)*log(2*pi)+0.5*log(obj.det_sigma);
        end
        
        function vals = draw(obj,n_draws)
            assert(any(size(obj.mu,1)==[1,n_draws]) && any(size(obj.sigma,3)==[1,n_draws]),...
                        'Obj must either have single value for parameters or the same number as wish to be sampled');
            if ismatrix(obj.chol_sigma)
               vals = bsxfun(@plus,obj.mu,randn(n_draws,size(obj.mu,2))*obj.chol_sigma);
            else
               vals = bsxfun(@plus,obj.mu,...
                                squeeze(sum(bsxfun(@times,randn(n_draws,size(obj.mu,2)),...
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
                
        function c = cdf(obj,vals) %#ok<INUSD,STOUT>
            error('Only write mv_gaussian_class cdf if required as less trivial');
        end
        
        function log_c = log_cdf(obj,vals) %#ok<INUSD,STOUT>
            error('Only write mv_gaussian_class log_cdf if required as less trivial');
        end
        
    end
end