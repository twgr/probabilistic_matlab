classdef discrete_class < base_primitive
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
        
        function vals = draw(obj,n_draws)
            assert(any(size(obj.p,1)==[1,n_draws]),...
                'Obj must either have single value for parameters or the same number as wish to be sampled');
            vals = sum(bsxfun(@gt,rand(n_draws,1),obj.cum_p),2)+1;
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
        
        function log_p = log_pdf(obj,vals)
            log_p = log(obj.pdf(vals));
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