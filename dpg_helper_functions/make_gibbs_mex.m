% MAKE_GIBBS_MEX   Generate MEX-function gibbs_update_j_mex from gibbs_update_j.
% 
% Script generated from project 'gibbs_update_j.prj' on 26-Jan-2016.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.TargetLang = 'C++';
cfg.GenerateReport = true;

%% Define argument types for entry-point 'gibbs_update_j'.
ARGS = cell(1,1);
ARGS{1} = cell(4,1);
ARGS{1}{1} = coder.typeof(0,[320  1],[1 0]);
ARGS{1}{2} = coder.typeof(0);
ARGS{1}{3} = coder.typeof(0);
ARGS{1}{4} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg gibbs_update_j -args ARGS{1}
