function log_Z = log_marginal(samples)
%log_marginal
%
% log_Z = log_marginal(samples)
%
% Provides an estimate for the marginal likelihood based on the full set of
% samples.  Currently only valid for smc inference.
%
% Tom Rainforth 05/07/16

assert(strcmpi(samples.options.inference_type,'smc'),'Only smc provides an overall marginal likehood estimate!')
log_Zs = sample.other_outputs.log_Zs;
max_log_Zs = max(log_Zs(:));
Zs = exp(log_Zs(:)-max_log_Zs);
log_Z = log(mean(Zs))+max_log_Zs;