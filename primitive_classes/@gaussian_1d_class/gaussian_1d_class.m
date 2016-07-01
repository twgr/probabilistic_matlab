classdef gaussian_1d_class < base_primitive
    properties
        mu
        sigma
    end
    
    methods
        
        function obj = gaussian_1d_class(mu,sigma)
            obj.mu      = mu;
            obj.sigma   = sigma;
        end
        
        function vals = draw(obj,n_draws)
            assert(any(size(obj.mu,1)==[1,n_draws]) && any(size(obj.sigma,1)==[1,n_draws]),...
                'Obj must either have single value for parameters or the same number as wish to be sampled');
            vals = randn(n_draws,1)*obj.sigma+obj.mu;
        end
        
        function log_p = log_pdf(obj,vals)
            log_p = bsxfun(@minus,-0.5*log(2*pi*obj.sigma.^2),...
                                  bsxfun(@rdivide,bsxfun(@minus,vals,obj.mu),...
                                                  (sqrt(2))*obj.sigma).^2);
        end
        
        function p = pdf(obj,vals)
            p = normpdf(vals,obj.mu,obj.sigma);
        end
                
        function c = cdf(obj,vals)
            c = normcdf(vals,obj.mu,obj.sigma);
        end
        
        function log_c = log_cdf(obj,vals)
            log_c = log(obj.cdf(vals));
        end
        
    end
end