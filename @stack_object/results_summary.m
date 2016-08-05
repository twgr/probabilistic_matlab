function summary = results_summary(samples,b_include_convergence)
%summary = convergence_summary(samples)
%
% Calculates the summary fields below for all the
% numeric fields stored in samples.var.
%
% summary fields:
%   ESS
%   Mean
%   Std_dev
%   Skewness
%   Excess_kurtosis
%
% Inputs:
%   samples
%   b_include_convergence = If true (default), calculates full convergence 
%                           results, else only the final results
%
% Outputs:
%   summary = Structure with a field for con, options, other_outputs and
%             var. The first three are taken directly from samples, 
%             the last is a structure of summaries for each of the numeric
%             variables in samples.var
%
% Tom Rainforth 05/08/16

if ~exist('b_include_convergence','var') || isempty(b_include_convergence)
    b_include_convergence = true;
end

output_field_names = {'ESS','Mean','Std_dev','Skewness','Excess_kurtosis'};

summary.other_outputs = samples.other_outputs;
summary.options = samples.options;
summary.con = samples.con;

var_fields = fields(samples.var);
b_numeric = cellfun(@(x) isnumeric(samples.var.(x)), var_fields);
var_fields = var_fields(b_numeric);

for n=1:numel(var_fields)
    
    var_name = var_fields{n};
    
    if b_include_convergence
        ESS = ess_convergence(samples,var_name); %#ok<NASGU>    
        [Mean,Std_dev,Skewness,Excess_kurtosis] = empirical_moments_convergence(samples,1:4,var_name); %#ok<ASGLU>
    else
        ESS = ess(samples,var_name); %#ok<NASGU>    
        [Mean,Std_dev,Skewness,Excess_kurtosis] = empirical_moments(samples,1:4,var_name); %#ok<ASGLU>
    end
    
    for nof = 1:numel(output_field_names)
        eval(['summary.var.' var_name '.' output_field_names{nof} ' = ' output_field_names{nof} ';']);
    end
end