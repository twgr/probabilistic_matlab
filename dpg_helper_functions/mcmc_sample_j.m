function [j_sample, estimated_marginals] = mcmc_sample_j(log_Z,P,n_steps,n_burn_in)
 
points = NaN(n_steps,P);

points(1,:) = randperm(numel(log_Z),P);
i_left = 1:numel(log_Z);
i_left(points(1,:)) = [];
i_include_proposals = randi(numel(i_left),n_steps-1,1);
i_exclude_proposals = randi(P,n_steps-1,1);

old_log_weight = sum(log_Z(points(1,:)));

for n=2:n_steps
    propose_points = points(n-1,:);
    propose_points(i_exclude_proposals(n-1)) = i_left(i_include_proposals(n-1));
    new_log_weight = sum(log_Z(propose_points));
        
    b_accept = rand<min(1,exp(new_log_weight-old_log_weight));
    
    if b_accept
        points(n,:) = propose_points;
        i_left(i_include_proposals(n-1)) = points(n-1,i_exclude_proposals(n-1));
        old_log_weight = new_log_weight;
    else
        points(n,:) = points(n-1,:);
    end
end

j_sample = points(end,:);
points = points(n_burn_in:end,:);
estimated_marginals = accumarray(points(:),ones(numel(points),1),[numel(log_Z),1]);

estimated_marginals = estimated_marginals/sum(estimated_marginals);
j_sample = j_sample(:);

end

