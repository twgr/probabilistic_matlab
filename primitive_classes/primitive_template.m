classdef primitive_template
 % Minimum requirement for making a new distribution primitive.  More
 % advanced usage is permitted (see current primitives) but all new
 % primitives must contain a minimum of this template.  Additional
 % properites and methods can be defined, but the methods defined here
 % (with the exception of the constructor, cdf and log_cdf) must take in 
 % and return the same arguments as given in this template.
    
    properties
        property_1
        property_2
        % ...
    end
    
    methods
        
        function obj = primitive_template(property_1,property_2)
            % Constructor.  The design of this is flexible but the returned
            % object must contain all the information required to execute
            % the sample and observe statements.  Different rows of the
            % property will correspond to different instances and so all
            % operations must support an arbitrary size for this first
            % dimension unless the property is contant or each instance is
            % a matrix.  If each instance is a matrix, see the
            % mv_gaussian_class for guidance.
            obj.property_1 = property_1;
            obj.property_2 = property_2;
        end
        
        function vals = sample(obj)
            global sample_size;
            
            % This assert is common and ensures that properties are either
            % inline with the size of the sample_size or have first
            % dimension 1
            assert(any(size(obj.property_1,1)==[1,sample_size]) && any(size(obj.property_2,1)==[1,sample_size]),...
                'Obj must either have single value for parameters or the same number as wish to be sampled');
            
            vals = ; % Sampled values
        end
        
        function log_p = observe(obj,vals)
            log_p = ; % Log likelihood of values
        end
        
        function p = pdf(obj,vals)
            p = ; % Probabilitiy denisty function (=exp(obj.observe,vals))
        end
        
        function c = cdf(obj,vals)
            c = ; % Cumulative density function.  Not currently required but preferred.
        end
        
        function log_c = log_cdf(obj,vals)
            log_c = ; % % Log cumulative density function.  Not currently required but preferred.
        end
        
    end
end