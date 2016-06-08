classdef gaussian_1d_class
    properties
        mu
        sigma
    end
    
    methods
        
        function obj = gaussian_1d_class(mu,sigma)
            obj.mu      = mu;
            obj.sigma   = sigma;
        end
        
        function vals = sample(obj)
            global sample_size;
            
            assert(any(size(obj.mu,1)==[1,sample_size]) && any(size(obj.sigma,1)==[1,sample_size]),...
                'Obj must either have single value for parameters or the same number as wish to be sampled');
            vals = randn(sample_size,1)*obj.sigma+obj.mu;
        end
        
        function log_p = observe(obj,vals)
            log_p = log(obj.pdf(vals));
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