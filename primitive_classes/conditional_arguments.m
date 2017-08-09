function [X, other_output] = conditional_arguments(X,argument_function,calling_function,varargin)

global sample_size;

N = sample_size;

args = argument_function(X);

[unique_args,~,argument_id] = unique(args,'rows');

N_u = size(unique_args,1);

X_template = stack_object;
X_template.con = X.con;
X_sub = repmat(X_template,N_u,1);
other_outputs = cell(N_u,1);

n_this = NaN(N_u,1);

for n_s = 1:N_u
    n_this(n_s) = sum(argument_id==n_s);
end

variables = fields(X.var);
assignment_order = cell(N_u,1);

for n_s = 1:N_u
    assignment_order{n_s} = find(argument_id==n_s);
    for n_v = 1:numel(variables)
        X_sub(n_s).var.(variables{n_v}) = X.var.(variables{n_v})(argument_id==n_s,:);
    end
    sample_size = n_this(n_s);
    [X_sub(n_s), other_outputs{n_s}] = calling_function(X_sub(n_s),unique_args(n_s,:),varargin{:});
    if size(other_outputs{n_s},1)==1
        other_outputs{n_s} = repmat(other_outputs{n_s},sample_size,1);
    end
end

sample_size = N;

[X.var, other_output] = struct_array_to_single_struct([X_sub(:).var], assignment_order, other_outputs);

end