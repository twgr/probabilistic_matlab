classdef gaussian_1d_class < base_primitive
    properties
        mu
        sigma
        Z      % Optional input, allows truncation of the distribution in
        % the cumulative normal space.
    end
    
    methods
        
        function obj = gaussian_1d_class(mu,sigma,Z)
            obj.mu      = mu;
            obj.sigma   = sigma;
            if exist('Z','var')
                obj.Z = Z;
            end
        end
        
        function obj = set.Z(obj,Z)
            if size(Z,2)==1
                Z = Z';
            end
            obj.Z = Z;
        end
        
        function vals = draw(obj,n_draws)
            if isempty(obj.Z)
                assert(any(size(obj.mu,1)==[1,n_draws]) && any(size(obj.sigma,1)==[1,n_draws]),...
                    'Obj must either have single value for parameters or the same number as wish to be sampled');
                vals = randn(n_draws,1)*obj.sigma+obj.mu;
            else
                Zs = rand(n_draws,1).*(obj.Z(:,2)-obj.Z(:,1))+obj.Z(:,1);
                vals = norminv(Zs,obj.mu,obj.sigma);
            end
        end
        
        function log_p = log_pdf(obj,vals)
            log_p = bsxfun(@minus,-0.5*log(2*pi*obj.sigma.^2),...
                bsxfun(@rdivide,bsxfun(@minus,vals,obj.mu),...
                (sqrt(2))*obj.sigma).^2);
            if ~isempty(obj.Z)
                log_p = log_p-log(obj.Z(:,2)-obj.Z(:,1));
            end
        end
        
        function p = pdf(obj,vals)
            p = exp(obj.log_pdf(vals));
        end
        
        function c = cdf(obj,vals)
            c = normcdf(vals,obj.mu,obj.sigma);
        end
        
        function log_c = log_cdf(obj,vals)
            log_c = log(obj.cdf(vals));
        end
        
    end
end