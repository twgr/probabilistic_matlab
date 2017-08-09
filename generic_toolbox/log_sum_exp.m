function [log_Z,w] = log_sum_exp(log_ws)

if all(log_ws==-inf)
    log_Z = -inf;
    w = ones(size(log_ws));
else
    z_max = max(log_ws);
    w = exp(log_ws-z_max);
    log_Z = z_max+log(sum(w));
end