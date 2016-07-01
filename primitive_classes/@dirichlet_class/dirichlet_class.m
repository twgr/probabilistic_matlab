classdef dirichlet_class < base_primitive
    properties
        alphas
    end
    
    properties (Hidden=true, SetAccess=private)
        log_beta_alpha
    end
    
    methods
        
        function obj = dirichlet_class(alphas)
            obj.alphas = alphas;
        end
        
        function obj = set.alphas(obj,alphas)
            obj.alphas = alphas;
            obj.log_beta_alpha = bsxfun(@minus,sum(gammaln(alphas),2),gammaln(sum(alphas,2))); %#ok<MCSUP>
        end
        
        function vals = draw(obj,n_draws)
            assert(any(size(obj.alphas,1)==[1,n_draws]),...
                'Obj must either have single value for parameters or the same number as wish to be sampled');
            
            if size(obj.alphas,1)==1
               vals = gamrnd(repmat(obj.alphas,n_draws,1),1,n_draws,size(obj.alphas,2));
            else
               vals = gamrnd(obj.alphas,1,n_draws,size(obj.alphas,2));
            end
            
            vals = bsxfun(@rdivide,vals,sum(vals,2));
            
        end
                
        function p = pdf(obj,vals)
            p = exp(obj.log_pdf(vals));
        end
        
        function log_p = log_pdf(obj,vals)
            log_p = bsxfun(@minus,sum(bsxfun(@times,obj.alphas-1,log(vals)),2),...
                                  obj.log_beta_alpha);
        end
        
        function c = cdf(obj,vals) %#ok<INUSD,STOUT>
            error('Dirichlet cdf unavailible');
        end
        
        function log_c = log_cdf(obj,vals) %#ok<INUSD,STOUT>
            error('Dirichlet cdf unavailible');
        end
        
    end
end