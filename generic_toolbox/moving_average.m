function Xm = moving_average(X,ni)
% function Xm = moving_average(X,ni)
%
% Calculates the moving average Xm, of the arbitrary dimension array X
% using the block sizes for ni
%
% Tom Rainforth 08/06/16

if ~all(mod(ni,2))
    error('The moving average sizes must be odd numbers otherwise they cannot centre about the correct point');
end

if (numel(ni)~=ndims(X)) && ~iscolumn(X) && ~isrow(X)
    error('ni needed for each of the dimensions of X');
end
    

Xm = NaN(size(X));
id = cell(1,numel(ni));
idsMean = cell(1,numel(ni));
idsAll = cell(1,numel(ni));

for n=1:prod(ni)
    [id{:}] = ind2sub(ni,n);    
    for m=1:numel(ni)
        idsMean{m} = (id{m}+floor(ni(m)/2)):ni(m):(size(X,m)-floor(ni(m)/2));
        idsAll{m} = id{m}:size(X,m);
    end        
    Xm(idsMean{:}) = blockwise_mean(X(idsAll{:}),ni);
end

end