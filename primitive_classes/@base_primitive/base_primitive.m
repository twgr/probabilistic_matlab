classdef base_primitive
% A base class from which the other primitives inherit.  Thus all required
% interupts and placed in here.  Currently does not do much other than set
% the sample size but using global variables and expanding the sample and
% observe functions allows significant extra functionality to be
% incorporated such as proposal adaptation etc.

    properties
    
    end
        
    methods
        function var = sample(obj)
            global sample_size
            
            var = obj.draw(sample_size);
        end
        
        function var = observe(obj,vals)
            var = obj.log_pdf(vals);
        end
        
    end    
end