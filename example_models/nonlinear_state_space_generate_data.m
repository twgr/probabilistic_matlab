function [sampling_functions,observe_functions] = nonlinear_state_space_generate_data(model_inputs) 
%% CHANGE THE VALUES ASSIGNED TO THE VARIABLES IN THIS SECTION

observations = model_inputs.observations;

if ~isfield(model_inputs,'X_1_mu')
    X_1_mu = 0;
else
    X_1_mu = model_inputs.X_1_mu;
end

if ~isfield(model_inputs,'X_1_sigma')
    X_1_sigma = sqrt(5);
else
    X_1_sigma = model_inputs.X_1_sigma;
end

if ~isfield(model_inputs,'Vn_mu')
    Vn_mu = 0;
else
    Vn_mu = model_inputs.Vn_mu;
end

if ~isfield(model_inputs,'Vn_sigma')
    Vn_sigma = sqrt(10);
else
    Vn_sigma = model_inputs.Vn_sigma;
end

if ~isfield(model_inputs,'Wn_sigma')
    Wn_sigma = sqrt(10);
else
    Wn_sigma = model_inputs.Wn_sigma;
end

n_observes = size(observations,1); 

observe_funcs_to_call = ones(n_observes,1); 

sample_funcs_to_call = [1;2*ones(n_observes-1,1)]; 

observe_args = num2cell(observations,2); 

sample_args = cell(n_observes,5);
sample_args(1,:) = {X_1_mu,X_1_sigma,Vn_mu,Vn_sigma,Wn_sigma};
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

function X = sample_1(X,X_1_mu,X_1_sigma,Vn_mu,Vn_sigma,Wn_sigma,varargin)

X.con.initial_distribution = gaussian_1d_class(X_1_mu,X_1_sigma);
X.con.rand_transition = gaussian_1d_class(Vn_mu,Vn_sigma);
X.con.likelihood = gaussian_1d_class(0,Wn_sigma);

X.var.x = X.con.initial_distribution.sample;
X.var.Y = (X.var.x.^2)/20+X.con.likelihood.sample;

end

function X = sample_2(X,n,varargin)

X.var.x = [X.var.x, transition_x_mean(X.var.x(:,end),n)+X.con.rand_transition.sample];
X.var.Y = [X.var.Y, (X.var.x(:,end).^2)/20+X.con.likelihood.sample];

end

%% USE THIS SECTION TO DEFINE HELPER FUNCTIONS

function log_w = observe_1(X,observation)

log_w = zeros(size(X.var.x,1),1);

end

function x_n_plus_1 = transition_x_mean(x_n,n)
    x_n_plus_1 = x_n/2+25*(x_n./(1+x_n.^2))+8*cos(1.2*n);
end