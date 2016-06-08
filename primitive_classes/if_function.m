function [X, other_output] = if_function(X,bool_test,true_return,false_return,varargin)
%if_function Function for if statements in probabilistic matlab
%
% X = if_function(X,bool_test,true_return,false_return,N)
%
% Inputs:  X         = Of stack_object type giving current state
%          bool_test = Annoynmous function that takes X and varargin and
%                      returns a
%                      boolean vector for the different instances stored in
%                      X.var, noting that the fields of X.var are arrays
%                      or cell arrays where different rows represent
%                      different instances.
%         true_return = Annoymous function that takes X and varargin as
%                       inputs and returns an updated X using one or more
%                       sample statements.  The syntax is as per the
%                       sample_n functions in the model file, except that
%                       changes to X.con are not allowed as by definition
%                       these must depended on the value of bool_test and
%                       so are not constants.
%                       Also allowed is a name - scalar value pair such as
%                       {'v',6};
%        false_return = As true_return but returned in the cases were
%                       bool_test is false.
%            varargin = Additional arguments that will be passed to the
%                       true_return and false_return functions.
%
% Output: X         = Updated stack_object with new samples added
%      other_output = Allows additional outputs to be returned such as a
%                     log weight for conditional observes.  Has the same
%                     restrictions as variables in X.var except that when a
%                     cell array it is allowed to have multiple columns
%                     corresponding to different outputs.

global sample_size;

N = sample_size;

bTrue = bool_test(X);

X_true = stack_object;
X_true.con = X.con;
X_false = stack_object;
X_false.con = X.con;

variables = fields(X.var);

for n_v = 1:numel(variables)
    X_true.var.(variables{n_v}) = X.var.(variables{n_v})(bTrue,:);
    X_false.var.(variables{n_v}) = X.var.(variables{n_v})(~bTrue,:);
end

n_true = sum(bTrue);
n_false = size(bTrue,1)-n_true;

if n_true>0
    sample_size = n_true;
    if iscell(true_return)
        X_true.var.(true_return{1}) = true_return{2}*ones(n_true,1);
        other_output_true = [];
    else
        [X_true, other_output_true] = true_return(X_true,varargin);
        if size(other_output_true,1)==1
            other_output_true = repmat(other_output_true,sample_size,1);
        end
    end
    if n_false==0
        X = X_true;
        other_output = other_output_true;
        sample_size = N;
        return;
    end
end

if n_false>0
    sample_size = n_false;
    if iscell(false_return)
        X_false.var.(false_return{1}) = false_return{2}*ones(n_false,1);
        other_output_false = [];
    else
        [X_false, other_output_false] = false_return(X_false,varargin);
        if size(other_output_false,1)==1
            other_output_false = repmat(other_output_false,sample_size,1);
        end
    end
    if n_true==0
        X = X_false;
        other_output = other_output_false;
        sample_size = N;
        return;
    end
end

sample_size = N;
i_true = find(bTrue);
i_false = find(~bTrue);
i_assign = [i_true;i_false];

X_dud = stack_object;
X_dud.con = X.con;

[X, other_output] = compose_two_sample_objects(X_dud,X_true,X_false,i_assign,...
                            n_true,n_false,other_output_true,other_output_false);


end