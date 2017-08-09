classdef beta_class < base_primitive
    properties
        a
        b
    end
    
    properties (Hidden=true, SetAccess=private)
        % This is the log normalization constant
        log_B
    end
    
    methods
        
        function obj = beta_class(a,b)
            obj.a = a;
            obj.b = b;
        end
        
        function obj = set.a(obj,a)
            obj.a = a;
            if ~isempty(obj.b) %#ok<MCSUP>
                obj.log_B = gammaln(obj.a)+gammaln(obj.b)-gammaln(obj.a+obj.b); %#ok<MCSUP>
            end
        end
        
        function obj = set.b(obj,b)
            obj.b = b;
            if ~isempty(obj.a) %#ok<MCSUP>
                obj.log_B = gammaln(obj.a)+gammaln(obj.b)-gammaln(obj.a+obj.b); %#ok<MCSUP>
            end
        end
        
        function vals = draw(obj,n_draws)
            assert(any(size(obj.a,1)==[1,n_draws]) && any(size(obj.b,1)==[1,n_draws]),...
                'Obj must either have single value for parameters or the same number as wish to be sampled');
            vals = betarnd(obj.a,obj.b,n_draws,1);
        end
        
        function p = pdf(obj,vals)
            p = exp(obj.log_pdf(vals));
        end
              
        function log_p = log_pdf(obj,vals)
            log_p = bsxfun(@minus,bsxfun(@times,obj.a-1,log(vals))+bsxfun(@times,obj.b-1,log(1-vals)),...
                                  obj.log_B);
        end
                
        function c = cdf(obj,vals)
            c = betacdf(vals,obj.a,obj.b);
        end
        
        function log_c = log_cdf(obj,vals)
            log_c = log(obj.cdf(vals));
        end
        
    end
end