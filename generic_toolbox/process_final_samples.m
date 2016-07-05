function samples = process_final_samples(samples_array,b_compress,T)
%process_final_samples
%
% Collapses the array of stack_objects produced by the inference methods to
% a single final object with well formatted and consistent weightings etc.
%
% Inputs:
%   samples_array = Array of stack_objects (n_iter x M)
%   b_compress = Whether to compress the results if not already compressed
%   T = Number of weighting function steps
%
% Outputs:
%   samples = Single stack_object storing details of all samples
%
% Tom Rainforth 08/06/16

samples_array = samples_array(:);

samples = stack_object;
samples.con = samples_array(1).con;

vars = [samples_array(:).var];
samples.var = struct_array_to_single_struct(vars,[]);

rw_cell = cell(numel(samples_array),1);
[rw_cell{:}] = deal(samples_array(:).relative_particle_weights);
rw_cell = cellfun(@transpose,rw_cell,'UniformOutput',false);
samples.relative_particle_weights = [rw_cell{:}]';


if b_compress
    samples = compress_samples(samples,T);
elseif ~isempty(samples_array(1).sparse_variable_relative_weights)
    if ~isstruct(samples_array(1).sparse_variable_relative_weights)
        [i,j,weight_vals] = find([samples_array(:).sparse_variable_relative_weights]);
        [sX1,sX2] = size(samples_array(1).sparse_variable_relative_weights);
        i = i+floor((j-1)/sX2)*sX1;
        j = mod(j,sX2);
        j(j==0) = sX2;
        samples.sparse_variable_relative_weights = sparse(i,j,weight_vals,numel(samples_array)*sX1,sX2);
    else
        sparse_variable_relative_weights = [samples_array(:).sparse_variable_relative_weights];
        field_names = fields(samples_array(1).sparse_variable_relative_weights);
        for n_f = 1:numel(field_names)
            [i,j,weight_vals] = find([sparse_variable_relative_weights(:).(field_names{n_f})]);
            [sX1,sX2] = size(samples_array(1).sparse_variable_relative_weights.(field_names{n_f}));
            i = i+floor((j-1)/sX2)*sX1;
            j = mod(j,sX2);
            j(j==0) = sX2;
            samples.sparse_variable_relative_weights.(field_names{n_f}) = sparse(i,j,weight_vals,numel(samples_array)*sX1,sX2);
        end
    end
end

end