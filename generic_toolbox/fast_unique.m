function [x_u, i_u, i_val] = fast_unique(X,bSorted)
%fast_unique Slightly faster version of matlabs unique that avoid overhead
%from argument checking and allows variable to be declared as already
%sorted
%
% Tom Rainforth

if exist('bSorted','var') && bSorted
    is = (1:numel(X))';
else    
    [X, is] = sort(X);
end
    
bDiff = [true;diff(X)~=0];

x_u = X(bDiff);

if nargout>1
    i_u = is(bDiff);
end

if nargout>2
    i_val = cumsum(bDiff);
    i_val(is) = i_val;
end

end