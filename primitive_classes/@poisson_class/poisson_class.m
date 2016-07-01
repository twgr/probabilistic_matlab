classdef poisson_class < base_primitive
    properties
        lambda
    end
        
    methods
        
        function obj = poisson_class(lambda)
            obj.lambda = lambda;
        end
        
        function vals = draw(obj,n_draws)
            assert(any(size(obj.lambda,1)==[1,n_draws]),...
                'Obj must either have single value for parameters or the same number as wish to be sampled');
            vals = poissrnd(obj.lambda,n_draws,1);
        end
        
        function log_p = log_pdf(obj,vals)
            bInt = round(vals)==vals;
            bInt = repmat(bInt,ceil(size(obj.lambda,1)/size(bInt,1)),ceil(size(obj.lambda,2)/size(bInt,2)));
            log_k_factorial = poisson_class.calc_log_k_factorial(vals);
            
            log_p = bsxfun(@minus,bsxfun(@minus,bsxfun(@times,vals,log(obj.lambda)),obj.lambda),log_k_factorial);
            log_p(~bInt) = -inf;          
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
    
    methods (Static)
        function log_k_factorial = calc_log_k_factorial(k)
            log_k_factorial = zeros(size(k));
            b_pos = k>0;
            if ~any(b_pos)
                return
            end
            log_k_factorial(b_pos) = poisson_class.calc_log_k_factorial(k(b_pos)-1)+log(k(b_pos));
        end
    end
end