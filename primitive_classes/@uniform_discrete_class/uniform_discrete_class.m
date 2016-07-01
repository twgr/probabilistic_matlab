classdef uniform_discrete_class < base_primitive
    properties
        minimum % Can have numerous columns
        maximum % Can have numerous columns
    end
    
    methods
        
        function obj = uniform_discrete_class(minimum,maximum)
            obj.minimum = minimum;
            obj.maximum = maximum;
        end
        
        function vals = draw(obj,n_draws)            
            assert(any(size(obj.minimum,1)==[1,n_draws]) && any(size(obj.maximum,1)==[1,n_draws]),...
                'Obj must either have single value for parameters or the same number as wish to be sampled');
            vals = obj.minimum+randi((obj.maximum-obj.minimum),n_draws,1);
        end
        
        function log_p = log_pdf(obj,vals)
            log_p = log(obj.pdf(vals));
        end
        
        function p = pdf(obj,vals)
            p = bsxfun(@times,1./(1+obj.maximum-obj.minimum),...
                bsxfun(@gt,vals,obj.minmum-1) & bsxfun(@lt,vals,obj.maximum+1));
        end
        
        function c = cdf(obj,vals)
            c = max(0,min(1,bsxfun(@rdivide,bsxfun(@minus,vals+1,obj.minimum),1+obj.maximum-obj.minimum)));
        end
        
        function log_c = log_cdf(obj,vals)
            log_c = log(obj.cdf(vals));
        end
        
    end
end