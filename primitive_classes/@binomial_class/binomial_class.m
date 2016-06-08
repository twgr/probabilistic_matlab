classdef binomial_class
    properties
        p % Can have numerous columns
        n % Can have numerous columns
    end
    
    methods
        
        function obj = binomial_class(p,n)
            obj.p = p;
            obj.n = n;
        end
        
        function vals = sample(obj)
            global sample_size;
            
            assert(any(size(obj.p,1)==[1,sample_size]) && any(size(obj.n,3)==[1,sample_size]),...
                'Obj must either have single value for parameters or the same number as wish to be sampled');
            vals = binornd(obj.n,obj.p,sample_size,1);
        end
        
        function log_p = observe(obj,vals)
            log_p = log(obj.pdf(vals));
        end
        
        function p = pdf(obj,vals)
            p = binopdf(vals,obj.n,obj.p);
        end
        
        function c = cdf(obj,vals)
            c = binocdf(vals,obj.n,obj.p);
        end
        
        function log_c = log_cdf(obj,vals)
            log_c = log(obj.cdf(vals));
        end
        
    end
end