function [summary, samples] = convergence_add_to_result_summaries(summary,samples)
%[summary, samples] = process_results(summary,samples)
% summary fields converted to convergence fields
%   Mean
%   Std_dev
%   Skewness
%   Excess_kurtosis

output_field_names = {'Mean','Std_dev','Skewness','Excess_kurtosis'};

variable_names = fields(summary);

for n=1:numel(variable_names)
    
    if iscell(variable_names{n})
        var_name = variable_names{n}{1};
        r1 = variable_names{n}{2};
    else
        var_name = variable_names{n};
        r1 = 1;
    end
        
    [Mean,Std_dev,Skewness,Excess_kurtosis] = empirical_moment_convergence(samples,1:4,var_name); %#ok<ASGLU>
    
    for nof = 1:numel(output_field_names)
        for nr = 1:r1
            if r1==1
                new_id = '';
            else
                new_id = ['_' num2str(nr)];
            end
            eval(['summary.' var_name new_id '.' output_field_names{nof} ' = ' output_field_names{nof} '(:,nr:r1:end);']);
        end
    end
    
    if r1~=1 && nargout>1
        X = samples.var.(var_name);
        for nr = 1:r1
            samples.var.([var_name '_' num2str(nr)]) = X(:,nr:r1:end);
        end
        samples.var = rmfield(samples.var,var_name);
        
        if isstruct(samples.sparse_variable_relative_weights)
            X = samples.sparse_variable_relative_weights.(var_name);
        else
            X = samples.sparse_variable_relative_weights;
        end
        for nr = 1:r1
            samples.sparse_variable_relative_weights.([var_name '_' num2str(nr)]) = X(:,nr:r1:end);
        end
        if isfield(samples.sparse_variable_relative_weights,var_name);
            samples.sparse_variable_relative_weights = rmfield(samples.sparse_variable_relative_weights,var_name);
        end
    end
end