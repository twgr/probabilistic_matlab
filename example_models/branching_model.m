function [sampling_functions,observe_functions] = branching_model(~)
%% CHANGE THE VALUES ASSIGNED TO THE VARIABLES IN THIS SECTION

n_observes = 1;
observe_funcs_to_call = ones(n_observes,1);
sample_funcs_to_call = ones(n_observes,1);
observe_args = cell(n_observes,1);
sample_args = cell(n_observes,1);


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

function X = sample_1(X,varargin)

% Sample variable l
X.con.count_prior = poisson_class(4);
X.var.r = X.con.count_prior.sample;

% Sample variable r
X = if_function(X,...
    @(X,~) X.var.r>4,...    % Bool for if statement
    {'l',6},...             % Return when true (name and value)
    @fib_branch);           % Return when false (function that updates X

end

function log_w = observe_1(X,varargin)

% Observe that poisson(X.var.l) = 6;
lik = poisson_class(X.var.l);
log_w = lik.observe(6);

end


%% Other required functions can be defined here

function [X,empty_output] = fib_branch(X,varargin)
empty_output = [];

% Set value of l to fib(3*r)+poissrnd(4);  Note that need to setup fib
% funciton to return a vector output for a vector input.
X.var.l = fib_vec(3*X.var.r)+X.con.count_prior.sample;
end


function a = fib_vec(ns)

% Memoized and vectorized code for returning an element of the fibonacci
% sequence

persistent a_vec

if isempty(a_vec)
    a_vec = NaN(40,1);
    for n=1:numel(a_vec)
        a_vec(n) = fib(n-1);
    end
elseif max(ns)>(numel(a_vec)-1)
    a_vec = [a_vec;NaN(max(ns)-(numel(a_vec)-1),1)];
    for n=numel(a_vec):max(ns)
        a_vec(n) = fib(n-1);
    end
end

a = a_vec(ns+1);

end


function a = fib(n)

% Recursive scalar calculation of fibonancci sequence

a = 0;
b = 1;
m = 0;

while m<n
    a_temp = b;
    b = a+b;
    a = a_temp;
    m = m+1;
end
end
