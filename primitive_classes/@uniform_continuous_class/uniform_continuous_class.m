classdef uniform_continuous_class
    properties
        minimum % Can have numerous columns
        maximum % Can have numerous columns
    end
    
    methods
        
        function obj = uniform_continuous_class(minimum,maximum)
            obj.minimum = minimum;
            obj.maximum = maximum;
        end
        
        function vals = sample(obj)
            global sample_size;
            
            assert(any(size(obj.minimum,1)==[1,sample_size]) && any(size(obj.maximum,1)==[1,sample_size]),...
                'Obj must either have single value for parameters or the same number as wish to be sampled');
            vals = obj.minimum+rand(sample_size,1).*(obj.maximum-obj.minimum);
        end
        
        function log_p = observe(obj,vals)
            log_p = log(obj.pdf(vals));            
        end
        
        function p = pdf(obj,vals)
            p = bsxfun(@times,1./(obj.maximum-obj.minimum),...
                                ones(max(size(vals,1),size(obj.minimum,1)),size(obj.minimum,2)));
        end
                
        function c = cdf(obj,vals)
            c = max(0,min(1,bsxfun(@rdivide,bsxfun(@minus,vals,obj.minimum),obj.maximum-obj.minimum)));
        end
        
        function log_c = log_cdf(obj,vals)
            log_c = log(obj.cdf(vals));
        end
        
    end
end