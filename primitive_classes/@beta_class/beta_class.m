classdef beta_class
    properties
        a
        b
    end
    
    methods
        
        function obj = beta_class(a,b)
            obj.a = a;
            obj.b = b;
        end
        
        function vals = sample(obj)
            global sample_size;
            
             assert(any(size(obj.a,1)==[1,sample_size]) && any(size(obj.b,1)==[1,sample_size]),...
                'Obj must either have single value for parameters or the same number as wish to be sampled');
            vals = betarnd(obj.a,obj.b,sample_size,1);
        end
        
        function log_p = observe(obj,vals)
            log_p = log(obj.pdf(vals));            
        end
        
        function p = pdf(obj,vals)
            p = betapdf(vals,obj.a,obj.b);
        end
                
        function c = cdf(obj,vals)
            c = betacdf(vals,obj.a,obj.b);
        end
        
        function log_c = log_cdf(obj,vals)
            log_c = log(obj.cdf(vals));
        end
        
    end
end