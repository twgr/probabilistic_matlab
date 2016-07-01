classdef gamma_class < base_primitive
    properties
        k       % shape parameter
        theta   % scale parameter - note this is 1/rate which is used in Anglican
    end
    
    methods
        
        function obj = gamma_class(k,theta)
            obj.k = k;
            obj.theta = theta;
        end
        
        function vals = draw(obj,n_draws)            
            assert(any(size(obj.k,1)==[1,n_draws]) && any(size(obj.theta,1)==[1,n_draws]),...
                'Obj must either have single value for parameters or the same number as wish to be sampled');
            
            vals = gamrnd(obj.k,obj.theta,n_draws,1);
        end
        
        function log_p = log_pdf(obj,vals)
            log_p = bsxfun(@minus,bsxfun(@minus,bsxfun(@times,obj.k-1,log(vals)),...
                bsxfun(@rdivide,vals,obj.theta)),...
                gammaln(obj.k)+obj.k.*log(obj.theta));
        end
        
        function p = pdf(obj,vals)
            p = gampdf(vals,obj.k,obj.theta);
        end
        
        function c = cdf(obj,vals)
            c = gamcdf(vals,obj.k,obj.theta);
        end
        
        function log_c = log_cdf(obj,vals)
            log_c = log(obj.cdf(vals));
        end
        
    end
end