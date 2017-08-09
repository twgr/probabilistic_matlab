classdef student_t_class < base_primitive
    properties
        nu
        mu
        sigma
    end
    
    methods
        
        function obj = student_t_class(nu,mu,sigma)
            obj.nu      = nu;
            obj.mu      = mu;
            obj.sigma   = sigma;
        end
        
        function vals = draw(obj,n_draws)
            assert(any(size(obj.mu,1)==[1,n_draws]) && any(size(obj.sigma,1)==[1,n_draws]) && any(size(obj.nu,1)==[1,n_draws]),...
                'Obj must either have single value for parameters or the same number as wish to be sampled');
            vals = bsxfun(@plus,bsxfun(@times,trnd(obj.nu,n_draws,1),obj.sigma),obj.mu);
        end
        
        function log_p = log_pdf(obj,vals)
            log_p = sum(bsxfun(@minus,gammaln((obj.nu+1)/2)-gammaln(obj.nu/2),...
                              bsxfun(@plus,log(sqrt(pi*obj.nu).*obj.sigma),...
                                bsxfun(@times,(obj.nu+1)/2,log(1+bsxfun(@times,1./obj.nu,...
                                                                             (bsxfun(@rdivide,(bsxfun(@minus,vals,obj.mu)),...
                                                                                        obj.sigma).^2)))))),2);
        end
        
        function p = pdf(obj,vals)
            p = exp(obj.log_pdf(vals));
        end
                
        function c = cdf(obj,vals)
            c = error('Cdf not currently provided');
        end
        
        function log_c = log_cdf(obj,vals)
            log_c = error('Cdf not currently provided');
        end
        
    end
end