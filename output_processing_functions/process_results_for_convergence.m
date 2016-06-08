function [summary, samples] = process_results_for_convergence(samples,varargin)
%[summary, samples] = process_results(samples,varargin)
% Varargin is each a cell of form:
%   name
%   number of different variables stored between columns
% or a string giving the name
%
% summary fields:
%   ESS
%   Mean
%   Std_dev
%   Skewness
%   Excess_kurtosis
%   (ESS_first_iter)
% These are given within a sub field of the name of the variable in
% question followed by _(dim) where (dim) is the dimension.
%
% If n_first_iter is provided (the number of samples in the first iteration)
% then ESS_first_iter is also returned as a summary field
%
% Returned samples is as per the input but with fields replaced buy the
% different reshape

output_field_names = {'ESS','Mean','Std_dev','Skewness','Excess_kurtosis'};

summary.other_outpus = samples.other_outputs;
summary.options = samples.options;
summary.con = samples.con;

n_particles = samples.options.n_particles;
if ~strcmpi('smc',samples.options.inference_type) && ~strcmpi('msmc',samples.options.inference_type) 
    n_iter = samples.options.n_iter;
    if isfield(samples.options,'n_nodes')
        n_nodes = samples.options.n_nodes;
    else
        n_nodes = 1;
    end
else
    n_nodes = 1;
    if isfield(samples.options,'n_nodes')
        n_iter = samples.options.n_nodes.*samples.options.n_iter;
    else
        n_iter = samples.options.n_iter;
    end
    samples.options.n_nodes = n_nodes;
    samples.options.n_iter = n_iter;
end


for n=1:numel(varargin)
    
    if iscell(varargin{n})
        var_name = varargin{n}{1};
        r1 = varargin{n}{2};
    else
        var_name = varargin{n};
        r1 = 1;
    end
    
    i_samples_first = bsxfun(@plus,(1:n_particles)',(n_particles*n_iter)*(0:1:(n_nodes-1)));
    ESS_first_iter = ess(samples,var_name,i_samples_first(:));
    if ~any(strcmpi(output_field_names,'ESS_first_iter'))
        output_field_names = [output_field_names, {'ESS_first_iter'}];
    end
    
    ESS = ess(samples,var_name); %#ok<NASGU>
    
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