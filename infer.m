function samples = infer(model_file_name,model_inputs,inference_type,varargin) %#ok<INUSL>

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

global_option_names = {'n_particles','n_iter','b_compress'};
global_default_values = {100, 100, 'default'};
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
        local_option_names = {'b_sparse'};
        local_default_values = {false};
        local_option_values = process_options(local_option_names,local_default_values,varargin{:});
                
        [samples, log_Z_total, log_Zs] = smc(sampling_functions,weighting_functions,global_option_values{:},local_option_values{:});
        other_outputs.log_Z_total = log_Z_total;
        other_outputs.log_Zs = log_Zs;
        
    case 'pimh'
        local_option_names = {'b_sparse'};
        local_default_values = {false};
        local_option_values = process_options(local_option_names,local_default_values,varargin{:});
        
        [samples, acceptance_ratio, log_Zs, b_accept] = pimh(sampling_functions,weighting_functions,global_option_values{:},local_option_values{:});
        other_outputs.acceptance_ratio = acceptance_ratio;
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
        
        [samples, log_Zs, iKept] = a_pgibbs(sampling_functions,weighting_functions,global_option_values{:},local_option_values{:});       
        other_outputs.log_Zs = log_Zs;
        other_outputs.iKept = iKept;
        
    case 'msmc'
        local_option_names = {'n_nodes','b_parallel','b_sparse'};
        local_default_values = {32, true, false};
        local_option_values = process_options(local_option_names,local_default_values,varargin{:});
        
        [samples, log_Z_total, log_Zs] = msmc(sampling_functions,weighting_functions,global_option_values{:},local_option_values{:});  
        other_outputs.log_Z_total = log_Z_total;
        other_outputs.log_Zs = log_Zs;
        
    case 'mpimh'
        local_option_names = {'n_nodes','b_parallel','b_sparse'};
        local_default_values = {32, true, false};
        local_option_values = process_options(local_option_names,local_default_values,varargin{:});
        
        [samples, acceptance_ratio, log_Zs] = mpimh(sampling_functions,weighting_functions,global_option_values{:},local_option_values{:});  
        other_outputs.acceptance_ratio = acceptance_ratio;
        other_outputs.log_Zs = log_Zs;
        
    case 'mpgibbs'
        local_option_names = {'n_nodes','b_parallel','initialization'};
        local_default_values = {32, true, []};
        local_option_values = process_options(local_option_names,local_default_values,varargin{:});
        
        [samples, log_Zs] = mpgibbs(sampling_functions,weighting_functions,global_option_values{:},local_option_values{:}); 
        other_outputs.log_Zs = log_Zs;
        
    case 'm_alt_pmcmc'
        local_option_names = {'n_nodes','b_parallel','initialization'};
        local_default_values = {32, true, []};
        local_option_values = process_options(local_option_names,local_default_values,varargin{:});
        
        [samples, log_Zs, iKept] = m_alt_pmcmc(sampling_functions,weighting_functions,global_option_values{:},local_option_values{:});
        other_outputs.log_Zs = log_Zs;       
        other_outputs.iKept = iKept;
      
    case 'mhalfhalf_pmcmc'
        local_option_names = {'n_nodes','b_parallel','n_csmc','initialization'};
        local_default_values = {32, true, 16, []};
        local_option_values = process_options(local_option_names,local_default_values,varargin{:});
        
        [samples, acceptance_ratios, log_Zs] = mhalfhalf_pmcmc(sampling_functions,weighting_functions,global_option_values{:},local_option_values{:});
        other_outputs.log_Zs = log_Zs;
        other_outputs.acceptance_ratios = acceptance_ratios;
        
    case 'ipmcmc'
        
        local_option_names = {'M','b_parallel','P','b_Rao_Black','n_conditional_gibbs_cycles','initial_retained_particles'};        
        local_default_values = {32, true, 16, true, 1, []};
        local_option_values = process_options(local_option_names,local_default_values,varargin{:});
        
        [samples, switching_rate, log_Zs, node_weights, sampled_indices] = ipmcmc(sampling_functions,weighting_functions,global_option_values{:},local_option_values{:}); 
        other_outputs.switching_rate = switching_rate;
        other_outputs.log_Zs = log_Zs;
        other_outputs.node_weights = node_weights;
        other_outputs.sampled_indices = sampled_indices;
                
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

for n=1:numel(option_names);
    bEqual = strcmpi(varargin,option_names{n});
    if any(bEqual)
        option_values{n} = varargin{find(bEqual,1)+1};
    end
end
end