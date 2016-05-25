function [i_resample, Z_step, n_times_sampled] = resample_step(log_weights,n_samples)

z_max = max(log_weights);
w = exp(log_weights(:)-z_max);
sum_w = sum(w);

drawsForResample = rand(n_samples,1);
w = w/sum(w);
edges = min([0;cumsum(w)],1);
edges(end) = 1;
[n_times_sampled,i_resample] = histc(drawsForResample,edges);
%[i_resample,i_r] = sort(i_resample);
%n_times_sampled = n_times_sampled(i_r);
%i_resample = datasample(1:numel(log_weights),n_samples,'Weights',w,'Replace',true)';
Z_step = z_max+log(sum_w)-log(numel(w));

end
