classdef bernoulli_class < base_primitive
    properties
        p % can have multiple columns to simultaenously encode many bernoulli distributions for vectorization
    end
    
    methods
        
        function obj = bernoulli_class(p)
            obj.p = p;
        end
        
        function vals = draw(obj,n_draws)
            
            assert(any(size(obj.p,1)==[1,n_draws]),...
                'Obj must either have single value for parameters or the same number as wish to be sampled');
            vals = double(bsxfun(@lt,rand(n_draws,size(obj.p,2)),obj.p));
        end
        
        function p = pdf(obj,vals)
            if size(obj.p,1)==1
                p = zeros(size(vals));
                p_vals = repmat(obj.p,size(vals,1),1);
                one_minus_p_vals = 1-p_vals;
                p(vals==1) = p_vals(vals==1);
                p(vals==0) = one_minus_p_vals(vals==0);
            elseif size(vals,1)==1
                p = zeros(size(obj.p));
                vals_rep = repmat(vals,size(p,1),1);
                p(vals_rep==1) = obj.p(vals_rep==1);
                p(vals_rep==0) = 1-obj.p(vals_rep==0);
            else
                p = zeros(size(vals));
                p(vals==1) = obj.p(vals==1);
                p(vals==0) = 1-obj.p(vals==0);
            end
        end
        
        function log_p = log_pdf(obj,vals)
            log_p = log(obj.pdf(vals));
        end
        
        function c = cdf(obj,vals)
            if size(obj.p,1)==1
                c = zeros(size(vals));
                one_minus_p_vals = repmat(1-obj.p,size(vals,1),1);
                c(vals>=1) = 1;
                b_zero_to_1 = vals>=0 & vals<1;
                c(b_zero_to_1) = one_minus_p_vals(b_zero_to_1);
            elseif size(vals,1)==1
                c = zeros(size(obj.p));
                vals_rep = repmat(vals,size(c,1),1);
                c(vals_rep>=1) = 1;
                b_zero_to_1 = vals_rep>=0 & vals_rep<1;
                c(b_zero_to_1) = 1-obj.p(b_zero_to_1);
            else
                c = zeros(size(vals));
                c(vals>=1) = 1;
                b_zero_to_1 = vals>=0 & vals<1;
                c(b_zero_to_1) = 1-obj.p(b_zero_to_1);
            end
        end
        
        function log_c = log_cdf(obj,vals)
            log_c = log(normcdf(vals,obj.mu,obj.sigma));
        end
        
    end
end