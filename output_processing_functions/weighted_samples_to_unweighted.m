function samples = weighted_samples_to_unweighted(samples,bRecompress)

p_fields = fields(samples.var);
n_samples = size(samples.var.(p_fields{1}),1);


if isprop(samples,'sparse_variable_relative_weights') && ~isempty(samples.sparse_variable_relative_weights)
    % Need to care that we don't blow the stack
    if ~exist('bRecompress','var')
        bRecompress = false;
    end    
    variable_sizes = NaN(numel(p_fields),2);
    for n_f = 1:numel(p_fields)
        variable_sizes(n_f,:) = size(samples.var.(p_fields{n_f}));
    end
    size_of_full_array = n_samples*sum(variable_sizes(:,2))*8; % In GB
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

if isfield(samples,'relative_weights')
    % Back compatibilitiy
    i_resample = datasample(1:n_samples,n_samples,'Weights',samples.relative_weights,'Replace',true)';
else
    i_resample = datasample(1:n_samples,n_samples,'Weights',samples.relative_particle_weights,'Replace',true)';
end


for n_f = 1:numel(p_fields)
    samples.var.(p_fields{n_f}) = samples.var.(p_fields{n_f})(i_resample,:);
end

if isfield(samples,'relative_weights')
    samples.relative_weights = [];
else
    samples.relative_particle_weights = [];
    samples.sparse_variable_relative_weights = [];
end

if bRecompress
    if numel(p_fields)==1
        n_resample_steps = size(samples.var.(p_fields{1}),2);
    else
        n_resample_steps = -1;
    end    
    samples = compress_samples(samples,1/n_samples*ones(n_samples,1),n_resample_steps);
end