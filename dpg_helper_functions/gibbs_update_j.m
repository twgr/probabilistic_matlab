function [j_sample, estimated_marginals] = gibbs_update_j(log_Z,P,n_runs,n_burn_in)

Z = exp(log_Z(:)-max(log_Z));
beta_sum = zeros(size(Z));

j_samples = NaN(n_runs,P);
j_samples(1,:) = 1:P;

for nr = 1:n_runs
    
    for p = 1:P
        beta_m = Z;
        beta_m(j_samples(nr,1:P~=p)) = 0;
        beta_m = beta_m/sum(beta_m);        
        cum_beta = cumsum(beta_m);
        cum_beta(end) = 1;
        j_this = 1+sum(rand>cum_beta);
        j_samples(nr,p) = j_this;
        
        if nr==1
            beta_sum = beta_sum+beta_m;
        end
    end
    
end

j_sample = j_samples(end,:);

if n_runs == 1
    estimated_marginals = beta_sum/sum(beta_sum); % This guards against numerical instability
else    
    j_sample = j_sample(n_burn_in+1:end,:);
    estimated_marginals = accumarray(j_sample(:),ones(numel(points),1),[numel(log_Z),1]);
    estimated_marginals = estimated_marginals/sum(estimated_marginals);
end

end