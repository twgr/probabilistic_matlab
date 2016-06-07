classdef stack_object
    
    properties
        options
        con
        var
        relative_particle_weights
        sparse_variable_relative_weights
        other_outputs
    end
    
    properties (Hidden=true)
        sparse_history
    end
    
end