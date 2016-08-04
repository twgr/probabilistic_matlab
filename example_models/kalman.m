function [sampling_functions,observe_functions] = kalman(model_inputs) 

if ~isfield(model_inputs,'omega1')
    omega1 = 0.7*pi;
else
    omega1 = model_inputs.omega1;
end

if ~isfield(model_inputs,'omega2')
    omega2 = 0.3*pi;
else
    omega2 = model_inputs.omega2;
end

if ~isfield(model_inputs,'omega3')
    omega3 = 0.05*pi;
else
    omega3 = model_inputs.omega3;
end

if ~isfield(model_inputs,'q')
    q = 1;
else
    q = model_inputs.q;
end

if ~isfield(model_inputs,'r')
    r = 0.1;
else
    r = model_inputs.r;
end

if ~isfield(model_inputs,'C')
    C = error('need to set C so have acess to it');
else
    C = model_inputs.C;
end

if ~isfield(model_inputs,'initial_state')
    initial_state = [0 1 1];
else
    initial_state = model_inputs.initial_state;
end

if ~isfield(model_inputs,'transition_scale')
    transition_scale = 0.99; %Needs to be less than 1 to stop the model being unstable
else
    transition_scale = model_inputs.transition_scale;
end

if ~isfield(model_inputs,'transition_mat')
    r1 = [cos(omega1),sin(-omega1),0; sin(omega1), cos(omega1), 0; 0, 0, 1];
    r2 = [cos(omega2),0,sin(-omega2); 0, 1, 0; -sin(omega2), 0, cos(omega2)];
    r3 = [1 0 0; 0 cos(omega3) -sin(omega3); 0 sin(omega3) cos(omega3)];
    transition_mat = transition_scale*(r3*r2*r1)'; % Transpose as we are doing x*R rather than R*x
else   
    transition_mat = model_inputs.transition_mat;
end

n_observes = size(model_inputs.observations,1); 

observe_funcs_to_call = ones(n_observes,1); 

sample_funcs_to_call = [1;2*ones(n_observes-1,1)]; 

observe_args = [num2cell(model_inputs.observations,2),num2cell((1:n_observes)',2),num2cell(size(transition_mat,1)*ones(n_observes,1),2)];

sample_args = [{q,r,C,initial_state,transition_mat};...
               cell(n_observes-1,5)];


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

function X = sample_1(X,q,r,C,initial_state,transition_mat,varargin)
    X.con.trans_noise = mv_gaussian_class(zeros(1,size(transition_mat,2)), q*eye(size(transition_mat,2)));
    X.con.obs_mat = C;
    X.con.lik = mv_gaussian_class(zeros(1,size(C,2)),r*eye(size(C,2)));
    X.con.trans_mat = transition_mat;
    X.var.x = bsxfun(@plus,initial_state(:)'*X.con.trans_mat,X.con.trans_noise.sample);
end

function X = sample_2(X,varargin)
    X.var.x = [X.var.x, bsxfun(@plus,X.var.x(:,end-2:end)*X.con.trans_mat,X.con.trans_noise.sample)];
end

function log_w = observe_1(X,observations,n_state,n_x_dims,varargin)
    log_w = X.con.lik.observe(bsxfun(@minus,X.var.x(:,(n_x_dims*(n_state-1)+1):n_x_dims*n_state)*X.con.obs_mat,observations));
end

%% USE THIS SECTION TO DEFINE HELPER FUNCTIONS
