function particles = compress_samples(particles, T)
%compress_samples
%
% Exploits the degeneracy caused by resampling to store the output using 
% sparse matrices and an implicit ancestral lineage.  This lineage is coded
% by the ordering of the samples - when a sample value is empty (signified
% by equalling zero in a the sparse array), it takes the value of the
% sample above (i.e. (i,j) takes its value from (i-1,j)).  If the sample
% above is also empty, it is equal to the sample above that and so on.
% This means that only samples that are unique at the current time step (or
% for discrete variables some previous time point) need to be stored,
% giving massive memory gains for long state sequences.  The particle
% weights are now also stored as a sparse array, will weights collapsed
% onto the sample where the variable values are provided.  The provided 
% output processing functions are overloaded to deal with this compressed
% format.  Care should be taken as a consequence of the coding, things such
% as naively taking the mean will give incorrect answers - the provided
% processing functions should be used instead or the compressed format
% avoided.
%
% Inputs:
%   particles = Uncompressed stack_object
%   T = Total number of weighting functions.  Note if the individual x_t
%       are multi-dimensional, this may be different to the array width
%
% Outputs:
%   particles = Compressed stack_object
%
% Tom Rainforth 08/06/16

p_fields = fields(particles.var);

b_zero_weight = particles.relative_particle_weights==0;
if all(b_zero_weight)
    particles.relative_particle_weights = 1e-100*ones(size(particles.relative_particle_weights));
else
    min_non_zero_weight = min(particles.relative_particle_weights(~b_zero_weight));
    weight_to_set_zeros_to = min(1e-100,exp(log(min_non_zero_weight)-32));
    if isnan(weight_to_set_zeros_to) || isinf(weight_to_set_zeros_to)
        weight_to_set_zeros_to = 1e-100;
    end
    particles.relative_particle_weights(b_zero_weight) = weight_to_set_zeros_to;
end

variable_sizes = NaN(numel(p_fields),2);
bContinuous = true(numel(p_fields),1);
bFullWidth = true(numel(p_fields),1);
for n_f = 1:numel(p_fields)
    assert(isnumeric(particles.var.(p_fields{n_f})),'Currently only supports compression when all variables are numeric');
    variable_sizes(n_f,:) = size(particles.var.(p_fields{n_f}));
    bContinuous(n_f) = rem(particles.var.(p_fields{n_f})(1),1)~=0;
    bFullWidth(n_f) = variable_sizes(n_f,2)==T;
end

iDiscrete = find(~bContinuous);

if T==1 && isempty(iDiscrete)
    warning('Compression does nothing if only single sample-weight step and continue variables');
    particles.relative_particle_weights = particles.relative_particle_weights;
    return
end

if (all(bFullWidth) && any(bContinuous)) || (numel(bFullWidth)==1)
    var_to_use = find(bContinuous);
    if isempty(var_to_use)
        var_to_use = 1;
    end
    [particles.var.(p_fields{var_to_use}),i_row_sort] = sortrows(particles.var.(p_fields{var_to_use}));
    w = particles.relative_particle_weights(i_row_sort);
    particles.relative_particle_weights = w;
    [sX1,sX2] = size(particles.var.(p_fields{var_to_use}));
    bUnique = [true(1,sX2);diff(particles.var.(p_fields{var_to_use}),[],1)~=0];
    if ~isempty(iDiscrete)
        % In the discrete case more care is required   
        for n_this = 2:sX2
            bUnique(:,n_this) = bUnique(:,n_this) | bUnique(:,n_this-1);
        end
    end
    indsU = find(bUnique);
    [iU,jU] = ind2sub([sX1,sX2],indsU);
    for n_f = 1:numel(p_fields)
        vars_this = particles.var.(p_fields{n_f})(indsU);
        b_zero = vars_this==0;
        if any(b_zero)
            warning('Zero values have been changed to 1e-32 in the compression');
            vars_this(b_zero) = 1e-32;
        end       
        particles.var.(p_fields{n_f}) = sparse(iU,jU,vars_this,sX1,sX2);
    end
    try
        cumsum_relative_weights = cumsum(w,'reverse');
    catch
        % In 2014 this format not supported
        wInv = w(end:-1:1);
        cumsum_relative_weights = cumsum(wInv);
        cumsum_relative_weights = cumsum_relative_weights(end:-1:1);
    end
    i_iU1 = iU==1;
    rel_w_1 = cumsum_relative_weights(iU);
    rel_w_2 = [cumsum_relative_weights(iU(2:end));0];
    rel_w_2([i_iU1(2:end);false]) = 0;
    weights_this = rel_w_1-rel_w_2;
    b_zero = weights_this==0;
    if any(b_zero)
        weights_this(b_zero) = min(1e-32,min(weights_this(~b_zero))/1e-10);
    end
    w = sparse(iU,jU,weights_this,sX1,sX2);
    particles.sparse_variable_relative_weights = w;
else
    % Here we are going to concatinate all the variables to do or
    % calculation to make sure variables are unique and then go back to
    % assign variables to the correct fields.  Care is taken to not store
    % the variables twice at any one point
    all_variables = NaN(variable_sizes(1,1),0);
    for n_f = 1:numel(p_fields)
        all_variables = [all_variables,particles.var.(p_fields{n_f})]; %#ok<AGROW>
        particles.var.(p_fields{n_f}) = [];
    end
    [all_variables,i_row_sort] = sortrows(all_variables);
    w = particles.relative_particle_weights(i_row_sort);
    particles.relative_particle_weights = w;
    [sX1,sX2] = size(all_variables);
    bUnique = [true(1,sX2);diff(all_variables,[],1)~=0];
    if ~isempty(iDiscrete)
        % In the discrete case more care is required   
        for n_this = 2:sX2
            bUnique(:,n_this) = bUnique(:,n_this) | bUnique(:,n_this-1);
        end
    end
    indsU = find(bUnique);
    [iU,jU] = ind2sub([sX1,sX2],indsU);
    try
        cumsum_relative_weights = cumsum(w,'reverse');
    catch
        % In 2014 this format not supported
        wInv = w(end:-1:1);
        cumsum_relative_weights = cumsum(wInv);
        cumsum_relative_weights = cumsum_relative_weights(end:-1:1);
    end
    for n_f = 1:numel(p_fields)
        bThis = jU<=variable_sizes(n_f,2);
        vars_this = all_variables(indsU(bThis));
        b_zero = vars_this==0;
        if any(b_zero)
            warning('Zero values have been changed to 1e-32 in the compression');
            vars_this(b_zero) = 1e-32;
        end
        particles.var.(p_fields{n_f}) = sparse(iU(bThis),jU(bThis),vars_this,sX1,variable_sizes(n_f,2));
        i_iU1 = iU(bThis)==1;
        rel_w_1 = cumsum_relative_weights(iU(bThis));
        rel_w_2 = [cumsum_relative_weights(iU([false;bThis(2:end)]));0];
        rel_w_2([i_iU1(2:end);false]) = 0;
        weights_this = rel_w_1-rel_w_2;
        b_zero = weights_this==0;
        if any(b_zero)
            weights_this(b_zero) = min(1e-32,min(weights_this(~b_zero))/1e-10);
        end
        particles.sparse_variable_relative_weights.(p_fields{n_f}) = sparse(iU(bThis),jU(bThis),weights_this,sX1,variable_sizes(n_f,2));
        
        iU = iU(~bThis);
        jU = jU(~bThis)-variable_sizes(n_f,2); 
        all_variables = all_variables(:,variable_sizes(n_f,2)+1:end);
        if n_f~=numel(p_fields)
            indsU = sub2ind([sX1,size(all_variables,2)],iU,jU);
        end
    end
end

end