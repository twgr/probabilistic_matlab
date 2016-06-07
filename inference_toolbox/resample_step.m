function [i_resample, Z_step, n_times_sampled] = resample_step(log_weights,n_samples)
%resample_step Samples with replacement from a set of log weights
%
% Inputs: 
%   log_weights = Log weights to resample from
%   n_samples = Number of samples to generate
%
% Outputs:
%   i_resample = Sampled indices
%   Z_step = sum(exp(log_weights))
%   n_times_sampled = Number of times each particle was resampled
%
% Tom Rainforth 07/06/16

z_max = max(log_weights);
w = exp(log_weights(:)-z_max);
sum_w = sum(w);

drawsForResample = rand(n_samples,1);
w = w/sum(w);
edges = min([0;cumsum(w)],1);
edges(end) = 1;
[n_times_sampled,i_resample] = histc(drawsForResample,edges);
Z_step = z_max+log(sum_w)-log(numel(w));

end
