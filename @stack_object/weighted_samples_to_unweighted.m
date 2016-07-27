function samples = weighted_samples_to_unweighted(samples,resample_method,n_samples_take,bRecompress)
%weighted_samples_to_unweighted
%
% samples = ...
%       weighted_samples_to_unweighted(samples,resample_method,bRecompress)
%
% Converts a stack_object with weighted samples to one where all the
% samples are unweighted by sampling with replacement using the resampling
% method defined by resample_method.
%
% Inputs
%   samples = Object of type stack_object
%   resample_method = Method for resampling, see resample_indices.m or
%             leave blank for default behaviour (stratified resampling)
%   n_samples_take = Number of samples to resample down to in the 
%             resampling. 
%                       Default = same number of samples as currently have
%   bRecompress = Whether to call compress_samples on the final samples.
%                       Default=false
%
% Outputs
%   samples = Object of type stack_object where all samples are equally
%             weighted.
%
% Tom Rainforth 05/07/16


if ~exist('resample_method','var')
    resample_method = []; % Take the default defined further on
end

p_fields = fields(samples.var);
n_samples_total = size(samples.var.(p_fields{1}),1);

if ~isempty(samples.sparse_variable_relative_weights)
    % Need to care that we don't blow the stack
    if ~exist('bRecompress','var') || isempty(bRecompress)
        bRecompress = false;
    end    
    variable_sizes = NaN(numel(p_fields),2);
    for n_f = 1:numel(p_fields)
        variable_sizes(n_f,:) = size(samples.var.(p_fields{n_f}));
    end
    size_of_full_array = n_samples_total*sum(variable_sizes(:,2))*8; % In GB
    try
        memory_stats = memory;
        largest_array = memory_stats.MaxPossibleArrayBytes;
    catch
        % memory function is only availible in windows
        largest_array = 4e9;
    end
    if size_of_full_array>(largest_array/3)
        error(sprintf('Converting to unweighted samples currently requires going through \nthe full sized array which would blow you RAM'));
    end
    for n_f = 1:numel(p_fields)
        samples.var.(p_fields{n_f}) = convert_to_full_array(samples.var.(p_fields{n_f}));
    end
else
    bRecompress = false;
end

if ~exist('n_samples_take','var') || isempty(n_samples_take)
    n_samples_take = n_samples_total;
end

samples.sparse_variable_relative_weights = [];
samples = resample_particles(samples, samples.relative_weights, n_samples_take, resample_method);

if bRecompress
    if numel(p_fields)==1
        n_resample_steps = size(samples.var.(p_fields{1}),2);
    else
        n_resample_steps = -1;
    end    
    samples = compress_samples(samples,n_resample_steps);
end