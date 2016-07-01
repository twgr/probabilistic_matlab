function samples = infer(model_file_name,model_inputs,inference_type,varargin) %#ok<INUSL>
%infer
%
% Entry function for performing inference using toolbox.  Operates on
% problems written in the required format by providing the path to the
% file.  Calling this model file generates the sampling_functions and
% weighting_functions cell array of anonymous functions that specify the
% model in the format required by the inference algorithms.
%
% Inputs:
%   model_file_name = Path to target model.  See example_run.m or
%                     example_models folder
%   model_inputs = Input as required by the model function.  If multiple
%                  inputs are required, these should be combined into a
%                  structure both here and in the model code
%   inference_type = Algorithm to use.  Currently supported:
%              'ipmcmc','smc','pimh','pgibbs','a_gibbs','indepedent_nodes'
%   varargin = Series of option name add values
%       - 'no_path_add' somewhere in the options means the path adding is
%          skipped
%       - Others should all be in pairs of name - value, for example 
%               'n_particles',1e5,'n_iter',1000
%       - Options for all algorithms, Allowed values (Default)
%              n_particles, +ve integere >=2 (1000)
%              resample_method, 'stratified' | 'multinomial' | 'systematic'
%                               | 'residual' ('stratified')
%              n_iter, +ve integere >=1 (100)
%              b_compress, 0 | false | 1 | true | 'default' | 'end_only'
%                          ('default', this looks at the other inputs and 
%                           estimates whether we should compress)
%       - For algorithm specific options see functions for individual
%         algorithms
%
% Outputs:
%   samples = Object of type stack_object that stores details of outputs
%             and forms input to most of postprocesssing functions
% 
% Tom Rainforth 20/06/16

i_path_add = find(strcmpi(varargin,'no_path_add'));
if isempty(i_path_add)
    addpath(genpath(regexprep(mfilename('fullpath'),'infer','')));
else
    i_input = 1:numel(varargin);
    i_input(i_path_add) = [];
    varargin = varargin(i_input);
end

other_outputs = struct;

eval(['[sampling_functions,weighting_functions] = ' model_file_name '(model_inputs);']);

global_option_names = {'n_particles','resample_method','n_iter','b_compress'};
global_default_values = {100, 'stratified', 100, 'default'};
global_option_values = process_options(global_option_names,global_default_values,varargin{:});
n_particles = global_option_values{strcmpi(global_option_names,'n_particles')};
n_iter = global_option_values{strcmpi(global_option_names,'n_iter')};
n_samples_to_generate = n_particles*n_iter*numel(sampling_functions);
n_b_compress = find(strcmpi(global_option_names,'b_compress'));

if n_samples_to_generate>2e7 && (numel(sampling_functions)~=1)
    
    try        
        memory_stats = memory;
        largest_array = memory_stats.MaxPossibleArrayBytes;        
    catch        
        % memory function is only availible in windows
        largest_array = 4e9;        
    end
        
    if (n_samples_to_generate*8)>(largest_array/50) && (ischar(global_option_values{n_b_compress}) || ~(global_option_values{n_b_compress}))
        if strcmpi(global_option_values{n_b_compress},'end_only') || (~ischar(global_option_values{n_b_compress}) && global_option_values{n_b_compress})
            warning('Turning on compression as at risk of causing a crash');
        end
        if ~strcmpi(global_option_values{n_b_compress},'force_off')
            global_option_values{n_b_compress} = true;
        end
    elseif (n_samples_to_generate*8)>(largest_array/500) && strcmpi(global_option_values{n_b_compress},'default')
        global_option_values{n_b_compress} = true;
    elseif (n_samples_to_generate*8)>(largest_array/2000) && strcmpi(global_option_values{n_b_compress},'default')
        global_option_values{n_b_compress} = 'end_only';
    end
    
end

if strcmpi(global_option_values{n_b_compress},'default')
    global_option_values{n_b_compress} = false;
end

if strcmpi(global_option_values{n_b_compress},'end_only')
    b_compress_end = true;
    global_option_values{n_b_compress} = false;
else
    b_compress_end = false;
end


switch inference_type
    case 'smc'
        local_option_names = {'b_parallel','n_islands','prop_sub_sample'};
        local_default_values = {true,1,1};
        local_option_values = process_options(local_option_names,local_default_values,varargin{:});
                
        [samples, log_Zs] = smc(sampling_functions,weighting_functions,global_option_values{:},local_option_values{:});
        other_outputs.log_Zs = log_Zs;
        
    case 'pimh'
        local_option_names = {'b_Rao_Black'};
        local_default_values = {true};
        local_option_values = process_options(local_option_names,local_default_values,varargin{:});
        
        [samples, log_Zs, b_accept] = pimh(sampling_functions,weighting_functions,global_option_values{:},local_option_values{:});
        other_outputs.log_Zs = log_Zs;
        other_outputs.b_accept = b_accept;
        
    case 'pgibbs'
        local_option_names = {'b_Rao_Black','initial_retained_particle'};
        local_default_values = {true,[]};
        local_option_values = process_options(local_option_names,local_default_values,varargin{:});
        
        [samples, log_Zs] = pgibbs(sampling_functions,weighting_functions,global_option_values{:},local_option_values{:});       
        other_outputs.log_Zs = log_Zs;
        
    case 'a_pgibbs'
        local_option_names = {'b_Rao_Black','initial_retained_particle'};
        local_default_values = {true,[]};
        local_option_values = process_options(local_option_names,local_default_values,varargin{:});
        
        [samples, log_Zs, b_accept] = a_pgibbs(sampling_functions,weighting_functions,global_option_values{:},local_option_values{:});       
        other_outputs.log_Zs = log_Zs;
        other_outputs.b_accept = b_accept;
        
    case 'indepedent_nodes'
        local_option_names = {'b_Rao_Black','b_parallel','Ms','initial_retained_particles'};
        local_default_values = {true, true, [32,0,0], []};
        local_option_values = process_options(local_option_names,local_default_values,varargin{:});
        
        [samples, log_Zs, b_accept] = mhalfhalf_pmcmc(sampling_functions,weighting_functions,global_option_values{:},local_option_values{:});
        other_outputs.log_Zs = log_Zs;
        other_outputs.b_accept = b_accept;
        
    case 'ipmcmc'
        
        local_option_names = {'b_Rao_Black','b_parallel','M','P','n_conditional_gibbs_cycles','initial_retained_particles'};        
        local_default_values = {true, true, 32, 16, 1, []};
        local_option_values = process_options(local_option_names,local_default_values,varargin{:});
        
        [samples, log_Zs, node_weights, sampled_indices, switching_rate] = ipmcmc(sampling_functions,weighting_functions,global_option_values{:},local_option_values{:}); 
        other_outputs.log_Zs = log_Zs;
        other_outputs.node_weights = node_weights;
        other_outputs.sampled_indices = sampled_indices;
        other_outputs.switching_rate = switching_rate;
                
    otherwise
        error('Invalid inference type.');
end

samples = process_final_samples(samples,b_compress_end,numel(weighting_functions));

options_pairs = [global_option_names, local_option_names, {'inference_type'};
                 global_option_values, local_option_values, {inference_type}];
samples.options = struct(options_pairs{:});
samples.other_outputs = other_outputs;
samples.options.inference_type = inference_type;

end

function option_values = process_options(option_names,default_values,varargin)
option_values = default_values;

% TODO add a check for spurious options 

% Set the options
for n=1:numel(option_names);
    bEqual = strcmpi(varargin,option_names{n});
    if any(bEqual)
        option_values{n} = varargin{find(bEqual,1)+1};
    end
end
end