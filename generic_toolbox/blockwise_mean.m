function Xm = blockwise_mean(X, ni)
% function Xm = blockwiseMean(X, ni)
%
% Calculates a blockwise mean of the arbitrary dimension array X using
% blocks of size given by the values of ni which must be a vector with the
% same number of values as the number of dimensions of X unless X is a 
% vector in which case ni can be a scalar. If X does not perfectly split
% into blocks of the correct size then the remainding data is simply
% ommited.
%
% e.g. Xm = blockwiseMean(randn(50,50,50),[10,3,2]);
%
% Tom Rainforth 08/06/16

if ndims(X)~=numel(ni)
    if iscolumn(X)
        ni = [ni,1];
    elseif iscolumn(X')
        ni = [1,ni];
    else
        error('Number of values in ni must equal the number of dimensions of X');
    end
end

ni = ni(:)';
S = size(X);
M = S - mod(S,ni);
Mn = M ./ ni;
reshapeInputArgs = [num2cell(ni); num2cell(Mn)];
XDims = arrayfun(@(x) 1:x, M, 'UniformOutput', false);
Xm = reshape(X(XDims{:}), reshapeInputArgs{:});
dimsToSum = -1+2*(1:numel(ni));
for n=1:numel(dimsToSum)
    Xm = sum(Xm, dimsToSum(n));
end
Xm = reshape(Xm / prod(ni), Mn);

end