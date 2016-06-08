function [sampling_functions,observe_functions] = model_template_commented(model_inputs) % CHANGE NAME TO FILE NAME
%% CHANGE THE VALUES ASSIGNED TO THE VARIABLES IN THIS SECTION
% This block done is by the user to define which functions are called and 
% their arguments.

% model_inputs is an argument that can be passed at the level of the infer
% function and allows for example data to be passed into the model.  It has
% not restrictions as its usage is completely controlled by the user

n_observes = ; % Change this to number of times an observe will be called.

observe_funcs_to_call = ; % List of observe function ids to be called.
% Must be a column vector of size n_observesx1.
% Note that as the function inputs can be varied between calls, the same
% observe function can be reused with different values for certain
% variables - for example the same likelihood function might be called
% numerous times with different data points.

sample_funcs_to_call = ; % Column vector of size n_observesx1 giving the
% number of the sample function to call prior to carrying out that observe.
% 0 indicates that nothing is sampled.  Note though that if nothing is
% sampled then it is more efficient to collapse the observe statements to a
% single statment as in this case resampling steps gain nothing.

observe_args = ; % Allows definition of arguments to the observe functions
% in addition to the stack of variables X. Must be a cell array with number 
% of rows equal to n_observes.  Columns correspond to the number of the
% input - e.g. for observe_1(X,a,b,varargin) then the value in the first
% column is assigned to variable a, the value of second column to b and so
% on.

sample_args = ;  % As observe_args but for the sample functions


%% LEAVE THIS SECTION AS IS

sampling_functions = cell(n_observes,1);
observe_functions = cell(n_observes,1);

for n=1:n_observes
    sampling_functions{n} = makefunction(true,sample_funcs_to_call(n),sample_args{n,:});
    observe_functions{n} = makefunction(false,observe_funcs_to_call(n),observe_args{n,:});
end

end

%% LEAVE THIS SECTION AS IS

function f = makefunction(bSample,number,varargin)

if bSample
    eval(['f = @(X) sample_' num2str(number) '(X,varargin{:});']);
else
    eval(['f = @(X) observe_' num2str(number) '(X,varargin{:});']);
end    
    
end

function X = sample_0(X,varargin) %#ok<DEFNU>
    % This is a placeholder for doing nothing;
end

%% WRITE THE SAMPLING FUNCTIONS AND OBSERVE FUNCTIONS IN THIS BLOCK

% Sample functions should all be named sample_n where n is a positive
% integer to idenify the function.  The values of sample_funcs_to_call will
% then correspond to the n that is called.  Probabilistic variables should
% be assigned to fields in X.var (for which there will be multiple
% instances) and constants to fields in X.con.  Note that X.con can be
% distribtuion objects provided they only have a single value of the
% parameter, e.g. X.con.p4 = poisson_class(4).
% 
% All functions should be of the form
%
% function X = sample_n(X,varargin) % n a positive integer
%   % Updates to X.  Note that the number of arguments passed is the width
%   % the sample_args cell array and therefore good practise is always to
%   % add varargin to the end of the function even if not used.
% end

% Observe functions should all be named observe_n where n is a positive
% integer to idenify the function.  The values of observe_funcs_to_call will
% then correspond to the n that is called.  All observe_n functions should
% be deterministic but branching is allowed using primitive/if_function.
% Note though that any updates to X will be temporary and will not persist
% outside the observe_n function.  Multiple observe calls to primitives are
% allowed, but the inference will treat these as a grouped single observe
% as resampling steps are only induced for every observe_n call and are based
% on their returned log_weight.
%
% All functions should be of the form
%
% function log_weight = observe_n(X,varargin) % n a positive integer
%   % Some temporary variable declarations.  Changes to X do not persist.
%   % The observe function of primitive classes can used to return the
%   % required log weight of the corresponding likelihood.  However, in
%   % theory any arbitrary likelihood function can be defined.
%   log_weight = ;
% end
%
% A small number of language constraints are:
%  1) Random draws should not be made other than with the supplied
%     primitives.  More complex distibituions must either be created
%     hierarchically using the provided primitives or by creating a new
%     primitive class.  The format for the primitive classes is consistent
%     so this should be easily done by adaptive a current primitive.
%  2) if and switch are allowed only for defining constants and must not
%     use sample, observe or be dependent on X.var.  Conditioning that
%     uses any of these should be done using the "primitives/if_function"
%  3) Functions using variables in X.var must be able to accept arguments
%     with an arbitrary length of first dimension (e.g. for numeric values
%     they must operate on both scalars and column vectors or on both row
%     vectors and matrices) and return outputs choose first dimension
%     corresponds to that of the input. For variables that whos dimension
%     may change (e.g. due to conditioning from "primitives/if_function"
%     will take the form of column cell arrays.  Again functions operating 
%     on such arbitrary dimension variables must operate on a Nx1 cell 
%     array where N is arbitrary. This condition should be easily 
%     satisified by either using vectorized code or by using a for loop 
%     around the first dimension of the variable in X.var and stacking the 
%     outputs.  Struct or object arrays are also allowed.



%% USE THIS SECTION TO DEFINE HELPER FUNCTIONS

% These have the same restirctions as the sample_n and observe_n functions.
% Although they may contain sample and observe statements, they must return
% and updated stack X if they do sample.