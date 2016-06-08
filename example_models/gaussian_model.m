function [sampling_functions,observe_functions] = gaussian_model(model_inputs)
%% CHANGE THE VALUES ASSIGNED TO THE VARIABLES IN THIS SECTION

if nargin==0 || isempty(model_inputs)
    model_inputs = [9;8]; % Default dataset if none provided
end

n_observes = 2; 

observe_funcs_to_call = ones(n_observes,1);         % Call the observe_1 function n_observes times
sample_funcs_to_call = [1;                          % First call sample_1
                        zeros(n_observes-1,1)];     % Sample nothing at later observe steps

observe_args = num2cell(model_inputs,2);            % model_inputs here is a column array of observed data points.
sample_args = [{1,  sqrt(5)};                       % Input 1 and sqrt(5) to first sample call
               repmat(cell(1,2),n_observes-1,1)];   % Input empties for any further calls                 

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

function X = sample_1(X,g_mean,g_std_dev,varargin)
    % Samples from a 1_d gaussian with mean g_mean and standard deviation
    % std_dev.  Assigns value to X.var.mu.  Note varargin must always be
    % added as other functions might take more arguments
    
    g1 = gaussian_1d_class(g_mean,g_std_dev);
    X.var.mu = g1.sample;
end


function log_w = observe_1(X,data_value)
	% Observe that normal(X.var.mu,sqrt(2))==data_value
    
    lik = gaussian_1d_class(X.var.mu, sqrt(2));
    log_w = lik.observe(data_value);
end


