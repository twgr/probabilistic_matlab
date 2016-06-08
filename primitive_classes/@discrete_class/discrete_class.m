classdef discrete_class
    properties
        p       % List of probabilities
    end
    
    properties (Hidden=true, SetAccess=private)
        cum_p   % Kept around to avoid the need for recalculation
    end
    
    methods
        
        function obj = discrete_class(p)
            obj.p = p;
        end
        
        function obj = set.p(obj,vals)
            obj.p = vals;
            obj.cum_p = cumsum(vals,2); %#ok<MCSUP>
        end
        
        function vals = sample(obj)
            global sample_size
            
            assert(any(size(obj.p,1)==[1,sample_size]),...
                'Obj must either have single value for parameters or the same number as wish to be sampled');
            vals = sum(bsxfun(@gt,rand(sample_size,1),obj.cum_p),2)+1;
        end
        
        function log_p = observe(obj,vals)
            log_p = log(obj.pdf(vals));
        end
        
        function p = pdf(obj,vals)
            if size(obj.p,1)==1
                p = obj.p(vals);
            elseif size(vals,1)==1
                p = obj.p(:,vals);
            else
                i = sub2ind(size(obj.p),(1:size(vals,1))',vals);
                p = obj.p(i);
                p = p(:);
            end
        end
        
        function c = cdf(obj,vals)
            if size(obj.cum_p,1)==1
                c = obj.cum_p(vals);
            elseif size(vals,1)==1
                c = obj.cum_p(:,vals);
            else
                i = sub2ind(size(obj.cum_p),(1:size(vals,1))',vals);
                c = obj.cum_p(i);
                c = c(:);
            end
        end
        
        function log_c = log_cdf(obj,vals)
            log_c = log(obj.cdf(vals));
        end
        
    end
end