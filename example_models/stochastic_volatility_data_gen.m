function [sampling_functions,observe_functions] = stochastic_volatility_data_gen(model_inputs) 
%% CHANGE THE VALUES ASSIGNED TO THE VARIABLES IN THIS SECTION

sigma = sqrt(0.9);
phi = 0.8;
beta = 0.7;
n_observes = 10000;

observe_funcs_to_call = ones(n_observes,1);
sample_funcs_to_call = [1;2*ones(n_observes-1,1)];
observe_args = cell(n_observes,1);

sample_args = cell(n_observes,3);
sample_args(1,:) = {sigma,phi,beta};
sample_args(2:end,1) = num2cell((2:n_observes)');


%% LEAVE THIS SECTION AS IS

sampling_functions = cell(n_observes,1);
observe_functions = cell(n_observes,1);

for n=1:n_observes
    sampling_functions{n} = makefunction(true,sample_funcs_to_call(n),sample_args{n,:});
    observe_functions{n} = makefunction(false,observe_funcs_to_call(n),observe_args{n,:});
end

end

function f = makefunction(bSample,number,varargin)

if bSample
    eval(['f = @(X) sample_' num2str(number) '(X,varargin{:});']);
else
    eval(['f = @(X) observe_' num2str(number) '(X,varargin{:});']);
end    
    
end

%% LEAVE THIS SECTION AS IS

function X = sample_0(X,varargin) %#ok<DEFNU>
    % This is a placeholder for doing nothing;
end

%% WRITE THE SAMPLING FUNCTIONS AND OBSERVE FUNCTIONS IN THIS BLOCK

function X = sample_1(X,sigma,phi,beta,varargin)

X.con.X_1_dist = gaussian_1d_class(0,sigma^2/(1-phi^2));
X.con.rand_transition = gaussian_1d_class(0,sigma);
X.con.W_dist = gaussian_1d_class(0,beta);
X.con.phi = phi;

X.var.x = X.con.X_1_dist.sample;
X.var.y = X.con.W_dist.sample*exp(X.var.x/2);

end

function X = sample_2(X,n_state,varargin)

X.var.x = [X.var.x, X.var.x(:,n_state-1)*X.con.phi+X.con.rand_transition.sample];
X.var.y = [X.var.y, X.con.W_dist.sample*exp(X.var.x(:,n_state)/2)];

end

%% USE THIS SECTION TO DEFINE HELPER FUNCTIONS

function log_w = observe_1(X,empty)

log_w = zeros(size(X.var.x,1),1);

end