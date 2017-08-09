function log_sum = log_minus_exp(log_1,log_2)

if log_1==-inf || log_2 == inf
    log_sum = -inf;
elseif log_1==inf || log_2 == -inf
    log_sum = inf;
else
    z_max = max([log_1,log_2]);
    w = exp([log_1,log_2]-z_max);
    log_sum = z_max+log(w(1)-w(2));
end