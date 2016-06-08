function [sampling_functions,observe_functions] = hmm(model_inputs) 
%% CHANGE THE VALUES ASSIGNED TO THE VARIABLES IN THIS SECTION

n_observes = size(model_inputs,1); 

observe_funcs_to_call = ones(n_observes,1); 

sample_funcs_to_call = [1;2*ones(n_observes-1,1)]; 

observe_args = [num2cell(model_inputs,2),num2cell((1:n_observes)'+1,2)];

sample_args = cell(n_observes,1);  


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

function X = sample_1(X,varargin)

X.con.initial_distribution = discrete_class(1/3*ones(1,3));
X.con.tr_1 = discrete_class([0.1 0.5 0.4]);
X.con.tr_2 = discrete_class([0.2 0.2 0.6]);
X.con.tr_3 = discrete_class([0.15 0.15 0.7]);
X.con.ob_1 = gaussian_1d_class(-1,1);
X.con.ob_2 = gaussian_1d_class(1,1);
X.con.ob_3 = gaussian_1d_class(0,1);

X.var.x = X.con.initial_distribution.sample;
X = sample_2(X,varargin);

end

function X = sample_2(X,varargin)

X = conditional_arguments(X,@(X) X.var.x(:,end), @transition_sample);

end

%% USE THIS SECTION TO DEFINE HELPER FUNCTIONS

function log_w = observe_1(X,observation,n_state)

[~,log_w] = conditional_arguments(X,@(X) X.var.x(:,n_state), @observe_state, observation);

end

function [X,other_outputs] = transition_sample(X,state,varargin)
    other_outputs = [];
    X.var.x = [X.var.x, X.con.(['tr_' num2str(state)]).sample];
end

function [X,log_w] = observe_state(X,state,data_point,varargin)
    log_w = X.con.(['ob_' num2str(state)]).observe(data_point);
end
