classdef poisson_class
    properties
        lambda
    end
    
    methods
        
        function obj = poisson_class(lambda)
            obj.lambda = lambda;
        end
        
        function vals = sample(obj)
            global sample_size;
            
            assert(any(size(obj.lambda,1)==[1,sample_size]),...
                'Obj must either have single value for parameters or the same number as wish to be sampled');
            vals = poissrnd(obj.lambda,sample_size,1);
        end
        
        function log_p = observe(obj,vals)
            log_p = log(obj.pdf(vals));            
        end
        
        function p = pdf(obj,vals)
            p = poisspdf(vals,obj.lambda);
        end
                
        function c = cdf(obj,vals)
            c = poisscdf(vals,obj.lambda);
        end
        
        function log_c = log_cdf(obj,vals)
            log_c = log(obj.cdf(vals));
        end
        
    end
end