function [summary, samples] = process_results(samples,n_particles,n_nodes,i_samples,varargin) %#ok<STOUT>
%[summary, samples] = process_results(samples,n_particles,n_nodes,i_samples,varargin)
% Varargin is each a cell of form:
%   name
%   number of different variables stored between columns
%
% summary fields:
%   ESS
%   Mean
%   Std_dev
%   Skewness
%   Excess_kurtosis
%   (ESS_first_iter)
% These are given within a sub field of the name of the variable in
% question followed by _(dim) where (dim) is the dimension.
%
% If n_first_iter is provided (the number of samples in the first iteration)
% then ESS_first_iter is also returned as a summary field
%
% Returned samples is as per the input but with fields replaced buy the
% different reshape

if iscell(n_particles)
    if ~exist('n_nodes','var')
        n_nodes = {};
    end
    if ~exist('i_samples','var')
        i_samples = {};
    end
    varargin = [{n_particles},{n_nodes},{i_samples},varargin];    
    n_particles = [];
    n_nodes = [];
    i_samples = [];
elseif iscell(n_nodes)
    if ~exist('i_samples','var')
        i_samples = {};
    end
    varargin = [{n_nodes},{i_samples},varargin];
    n_nodes = 1;
    i_samples = [];
elseif iscell(i_samples)
    varargin = [{i_samples},varargin];
    i_samples = [];
end

output_field_names = {'ESS','Mean','Std_dev','Skewness','Excess_kurtosis'};

for n=1:numel(varargin)
    
    if iscell(varargin{n})
        var_name = varargin{n}{1};
        r1 = varargin{n}{2};
    else
        var_name = varargin{n};
        r1 = 1;
    end
    
    if ~isempty(n_particles)
        n_samples = size(samples.var.(var_name),1);
        n_iter = n_samples/(n_nodes*n_particles);
        i_samples_first = bsxfun(@plus,(1:n_particles)',(n_particles*n_iter)*(0:1:(n_nodes-1)));
        ESS_first_iter = ess(samples,var_name,i_samples_first(:));
        if ~any(strcmpi(output_field_names,'ESS_first_iter'))
            output_field_names = [output_field_names, {'ESS_first_iter'}];
        end
    end
    
    ESS = ess(samples,var_name,i_samples); %#ok<NASGU>
    [Mean,Std_dev,Skewness,Excess_kurtosis] = empirical_moments(samples,1:4,i_samples,var_name); %#ok<ASGLU>
    
    for nof = 1:numel(output_field_names)
        eval([output_field_names{nof} ' = reshape(' output_field_names{nof} ',r1,[]);']);
        for nr = 1:r1
            if r1==1
                new_id = '';
            else
                new_id = ['_' num2str(nr)];
            end
            eval(['summary.' var_name new_id '.' output_field_names{nof} ' = ' output_field_names{nof} '(nr,:);']);
        end
    end
    if r1~=1 && nargout>1
        X = samples.var.(var_name);
        for nr = 1:r1
            samples.var.([var_name '_' num2str(nr)]) = X(:,nr:r1:end);
        end
        samples.var = rmfield(samples.var,var_name);
        
        if isstruct(samples.sparse_variable_relative_weights)
            X = samples.sparse_variable_relative_weights.(var_name);
        else
            X = samples.sparse_variable_relative_weights;
        end
        for nr = 1:r1
            samples.sparse_variable_relative_weights.([var_name '_' num2str(nr)]) = X(:,nr:r1:end);
        end
        if isfield(samples.sparse_variable_relative_weights,var_name);
            samples.sparse_variable_relative_weights = rmfield(samples.sparse_variable_relative_weights,var_name);
        end
    end
end