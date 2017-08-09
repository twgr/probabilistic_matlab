classdef uniform_continuous_class < base_primitive
    properties
        minimum % Can have numerous columns
        maximum % Can have numerous columns
        Z      % Optional input, allows truncation of the distribution in
        % the cumulative normal space.
    end
    
    methods
        
        function obj = uniform_continuous_class(minimum,maximum,Z)
            obj.minimum = minimum;
            obj.maximum = maximum;
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
            assert(any(size(obj.minimum,1)==[1,n_draws]) && any(size(obj.maximum,1)==[1,n_draws]),...
                'Obj must either have single value for parameters or the same number as wish to be sampled');
            if isempty(obj.Z)
                vals = obj.minimum+rand(n_draws,1).*(obj.maximum-obj.minimum);
            else
                Zs = obj.Z(:,1)+rand(n_draws,1).*(obj.Z(:,2)-obj.Z(:,1));
                vals = obj.minimum+Zs.*(obj.maximum-obj.minimum);
            end
        end
        
        function log_p = log_pdf(obj,vals)
            log_p = log(obj.pdf(vals));            
        end
        
        function p = pdf(obj,vals)
            p = bsxfun(@times,1./(obj.maximum-obj.minimum),...
                                ones(max(size(vals,1),size(obj.minimum,1)),size(obj.minimum,2)));
            if ~isempty(obj.Z)
                p = p./(obj.Z(:,2)-obj.Z(:,1));
            end
        end
                
        function c = cdf(obj,vals)
            c = max(0,min(1,bsxfun(@rdivide,bsxfun(@minus,vals,obj.minimum),obj.maximum-obj.minimum)));
        end
        
        function log_c = log_cdf(obj,vals)
            log_c = log(obj.cdf(vals));
        end
        
    end
end