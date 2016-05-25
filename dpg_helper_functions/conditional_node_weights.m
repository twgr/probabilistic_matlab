function w = conditional_node_weights(log_Z,P)
    
    Z = exp(log_Z(:)'-max(log_Z));
    
    [~,iMax] = max(Z);
    
    a = 1;
    for t=1:(P-1);
        a = recurse_a(a,Z);
    end
    a = a(:);
    
    if isnumeric(Z)
        w = sum(bsxfun(@times,a,bsxfun(@(x,y) x.^y,Z,(1:numel(a))')),1);
    else
        for n=1:numel(Z)
            w(n) = sum(a.*(Z(n).^reshape(1:numel(a),size(a))));
        end
    end
    
    if isnan(sum(w)) || sum(w)==0
        w = zeros(size(w));
        w(iMax) = 1;
    else
        w = w/sum(w);
    end
    
%     log_w2 = crazy_recursion(a,log_Z);
%     max_log_w2 = max(log_w2);
%     w2 = exp(log_w2-max_log_w2);
%     w2 = w2/sum(w2);
%     
%     disp(max(abs(w2-w)));
    
end