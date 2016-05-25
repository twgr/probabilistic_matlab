function [samples, relative_weights, switching_rate, log_Zs, node_weights, sampled_indices] = ...
    dpg(sampling_functions,weighting_functions,n_particles,...
        n_iter,b_compress,n_nodes,b_parrallel,n_csmc,...
        b_keep_samples_from_all_nodes,b_keep_all_samples_per_node,node_sampler,node_sampler_steps,node_sampler_burn_in,initialization)

if b_keep_samples_from_all_nodes
    if b_keep_all_samples_per_node
        relative_weights = NaN(n_particles*n_iter,n_nodes);
    else
        relative_weights = ones(n_iter,n_nodes);
    end
else
    if b_keep_all_samples_per_node
        relative_weights = NaN(n_particles*n_iter,1);
    else
        relative_weights = ones(n_iter,1);
    end
end

log_Zs = NaN(n_iter,n_nodes);
node_weights = NaN(n_iter,n_nodes);
sampled_indices = NaN(n_iter,n_csmc);

particles_iter = cell(n_nodes,1);
log_node_Z = NaN(n_nodes,1);
relative_weights_iter = cell(n_nodes,1);

if isempty(initialization)

if b_parrallel    
    parfor node = 1:n_nodes
        [particles_iter{node},log_node_Z(node),relative_weights_iter{node}] = smc_sweep(sampling_functions,weighting_functions,n_particles,false);
    end    
else    
    for node = 1:n_nodes
        [particles_iter{node},log_node_Z(node),relative_weights_iter{node}] = smc_sweep(sampling_functions,weighting_functions,n_particles,false);
    end    
end

else
    if iscell(initialization)
        initialization = [initialization{:}]';
    end
            
    if b_parrallel
        retained_particles_temp = [initialization;repmat(initialization(1),n_nodes-n_csmc,1)];
        parfor node = 1:n_nodes
            if node<=n_csmc
                [particles_iter{node},log_node_Z(node),relative_weights_iter{node}] = pg_sweep(sampling_functions,weighting_functions,n_particles,retained_particles_temp(node),false);
            else
                [particles_iter{node},log_node_Z(node),relative_weights_iter{node}] = smc_sweep(sampling_functions,weighting_functions,n_particles,false);
            end
        end
        clear retained_particles_temp
    else
        for node = 1:n_nodes
            if node<=n_csmc
                [particles_iter{node},log_node_Z(node),relative_weights_iter{node}] = pg_sweep(sampling_functions,weighting_functions,n_particles,initialization(node),false);
            else
                [particles_iter{node},log_node_Z(node),relative_weights_iter{node}] = smc_sweep(sampling_functions,weighting_functions,n_particles,false);
            end
        end
    end    
    
end

%[js, marginal_node_weights] = sample_conditional_nodes(log_node_weights,n_csmc);

[js, marginal_node_weights] = node_sampler(log_node_Z(:),n_csmc,node_sampler_steps,node_sampler_burn_in);

log_Zs(1,:) = log_node_Z;
node_weights(1,:) = marginal_node_weights;
sampled_indices(1,:) = js;

switching_rate = NaN(n_iter,1);
switching_rate(1) = 1;

ret_template = stack_object;
ret_template.con = particles_iter{1}.con;
retained_particles = repmat(ret_template,n_csmc,1);
p_fields = fields(particles_iter{1}.var);

for m=1:n_csmc
    
    ks = resample_step(max(log(relative_weights_iter{js(m)}),-1e4), 1);
    for n_f = 1:numel(p_fields)
        retained_particles(m).var.(p_fields{n_f}) = particles_iter{js(m)}.var.(p_fields{n_f})(ks,:);
    end
    
end

if ~b_keep_all_samples_per_node
    marginal_node_weights = marginal_node_weights/sum(marginal_node_weights(js));
end

returned_samples_ids = [js(:);setdiff((1:n_nodes)',js(:))];

for m=1:n_nodes
    
    if b_keep_samples_from_all_nodes
        if b_keep_all_samples_per_node
            samples(1,m) = particles_iter{returned_samples_ids(m)}; %#ok<AGROW>
            relative_weights(1:n_particles,m) = relative_weights_iter{returned_samples_ids(m)}*marginal_node_weights(returned_samples_ids(m));            
            if b_compress
                [samples(1,m), relative_weights(1:n_particles,m)] = compress_samples(samples(1,m),relative_weights(1:n_particles,m),numel(weighting_functions)); %#ok<AGROW>
            end
        elseif m<=n_csmc
            samples(1,m) = retained_particles(m); %#ok<AGROW>
            relative_weights(1,m) = marginal_node_weights(returned_samples_ids(m));
        else
            continue
        end
    elseif m==1
        if b_keep_all_samples_per_node
            samples(1,m) = particles_iter{returned_samples_ids(m)}; %#ok<AGROW>
            relative_weights(1:n_particles,m) = relative_weights_iter{returned_samples_ids(m)};
            if b_compress
                [samples(1,m), relative_weights(1:n_particles,m)] = compress_samples(samples(1,m),relative_weights(1:n_particles,m),numel(weighting_functions)); %#ok<AGROW>
            end
        else
            samples(1,m) = retained_particles(m); %#ok<AGROW>
            relative_weights(1,m) = 1;
        end
    else
        % Continue
    end
end

if ~b_compress && b_keep_all_samples_per_node
    S = whos('samples');
    s_mem = S.bytes*n_iter;
    if s_mem>5e7
        
        try
            memory_stats = memory;
            largest_array = memory_stats.MaxPossibleArrayBytes;
        catch
            % memory function is only availible in windows
            largest_array = 4e9;
        end
        
        if S.bytes*n_iter > (largest_array/20)
            warning('In danger of swamping memory and crashing, turning b_compress on');
            b_compress = true;
            for m=1:size(samples,2)
                [samples(1,m), relative_weights(1:n_particles,m)] = compress_samples(samples(1,m), relative_weights(1:n_particles,m), numel(sampling_functions));
            end
        end
    end
end

samples = repmat(samples,n_iter,1);

for iter=2:n_iter
    
    if b_parrallel
        retained_particles_temp = [retained_particles;repmat(retained_particles(1),n_nodes-n_csmc,1)];        
        parfor node = 1:n_nodes
            if node<=n_csmc
                [particles_iter{node},log_node_Z(node),relative_weights_iter{node}] = pg_sweep(sampling_functions,weighting_functions,n_particles,retained_particles_temp(node),false);
            else
                [particles_iter{node},log_node_Z(node),relative_weights_iter{node}] = smc_sweep(sampling_functions,weighting_functions,n_particles,false);
            end
        end
        clear retained_particles_temp
    else
        for node = 1:n_nodes
            if node<=n_csmc
                [particles_iter{node},log_node_Z(node),relative_weights_iter{node}] = pg_sweep(sampling_functions,weighting_functions,n_particles,retained_particles(node),false);
            else
                [particles_iter{node},log_node_Z(node),relative_weights_iter{node}] = smc_sweep(sampling_functions,weighting_functions,n_particles,false);
            end
        end
    end
        
    %[js, marginal_node_weights] = sample_conditional_nodes(log_node_weights,n_csmc);

    [js, marginal_node_weights] = node_sampler(log_node_Z(:),n_csmc,node_sampler_steps,node_sampler_burn_in);
        
    log_Zs(iter,:) = log_node_Z;
    node_weights(iter,:) = marginal_node_weights;    
    sampled_indices(iter,:) = js;
    switching_rate(iter) = sum(js>n_csmc)/n_csmc;
    
    ret_template = stack_object;
    ret_template.con = particles_iter{1}.con;
    retained_particles = repmat(ret_template,n_csmc,1);
    p_fields = fields(particles_iter{1}.var);
    
    for m=1:n_csmc
        ks = resample_step(max(log(relative_weights_iter{js(m)}),-1e4), 1);
        for n_f = 1:numel(p_fields)
            retained_particles(m).var.(p_fields{n_f}) = particles_iter{js(m)}.var.(p_fields{n_f})(ks,:);
        end
    end
    
    if ~b_keep_all_samples_per_node
        marginal_node_weights = marginal_node_weights/sum(marginal_node_weights(js));
    end
    
    returned_samples_ids = [js(:);setdiff((1:n_nodes)',js(:))];
    
    for m=1:n_nodes
        if b_keep_samples_from_all_nodes
            if b_keep_all_samples_per_node
                samples(iter,m) = particles_iter{returned_samples_ids(m)};
                relative_weights((1+(iter-1)*n_particles):(n_particles*iter),m) = relative_weights_iter{returned_samples_ids(m)}*marginal_node_weights(returned_samples_ids(m));
                if b_compress
                    [samples(iter,m), relative_weights((1+(iter-1)*n_particles):(n_particles*iter),m)]...
                                = compress_samples(samples(iter,m),relative_weights((1+(iter-1)*n_particles):(n_particles*iter),m),numel(weighting_functions));
                end
            elseif m<=n_csmc
                samples(iter,m) = retained_particles(m);
                relative_weights(iter,m) = marginal_node_weights(returned_samples_ids(m));
            else
                continue
            end
        elseif m==1
            if b_keep_all_samples_per_node
                samples(iter) = particles_iter{returned_samples_ids(m)};
                relative_weights((1+(iter-1)*n_particles):(n_particles*iter),m) = relative_weights_iter{returned_samples_ids(m)};
                if b_compress
                    [samples(iter,m), relative_weights((1+(iter-1)*n_particles):(n_particles*iter),m)]...
                        = compress_samples(samples(iter,m),relative_weights((1+(iter-1)*n_particles):(n_particles*iter),m),numel(weighting_functions));
                end
            else
                samples(iter) = retained_particles(m);
                relative_weights(iter,m) = 1;
            end
        else
            % Continue
        end
        
    end
    
end



end