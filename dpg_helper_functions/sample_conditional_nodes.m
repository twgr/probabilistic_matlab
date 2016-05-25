function [j_s,marginal_weights] = sample_conditional_nodes(log_Z,P)

j_s = NaN(P,1);
marginal_weights = NaN(numel(log_Z),1);

% First cream off the top and bottom 
[log_Z, iZOrder] = sort(log_Z,'descend');
i_unorder = sort(iZOrder);

log_Z_P = log_Z(P);
log_Z_P_plus_1 = log_Z(P);

i_def_in = find(log_Z_P>log_Z_P_plus_1+23); % If more than 10 orders of magnitude larger than the Pth plus 1 largest weight
i_def_not_in = find(log_Z_P<(log_Z_P-23));  % If less that 10 orders of magnitude less than the Pth largest weight

j_s(1:numel(i_def_in)) = i_unorder(i_def_in);
n_start = 1+numel(i_def_in);
marginal_weights(i_unorder(i_def_in)) = 1/P;
marginal_weights(i_unorder(i_def_not_in)) = 0;

i_left = setdiff(setdiff(1:numel(log_Z),i_def_in),i_def_not_in);
scale_left = 1-sum(marginal_weights);



for n=n_start:P
    
  w = conditional_node_weights(log_Z(i_left),P+1-n);
  if n==1
      marginal_weights(iZOrder(i_left)) = scale_left*w/sum(w);
  end
  cum_w = cumsum(w);
  j_local = i_left(1+sum(rand>cum_w));
  j_s(n) = i_unorder(j_local);
  i_left = setdiff(i_left,j_local);
  
end





