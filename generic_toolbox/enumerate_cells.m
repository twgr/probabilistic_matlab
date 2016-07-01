function E = enumerate_cells(varargin)
%enumerate_cells
%   E = enumerate_cells(varargin)
%
% Takes in an arbitrary number of inputs and enumerates into a single cell
% array with numel(varargin) dimensions.
%
% For example:
% % E = enumerate_cells([1,2],[4,5],[7,8])
% 
% E(:,:,1) = 
% 
%     [1x3 double]    [1x3 double]
%     [1x3 double]    [1x3 double]
% 
% 
% E(:,:,2) = 
% 
%     [1x3 double]    [1x3 double]
%     [1x3 double]    [1x3 double]
% 
% E{:}
% 
% ans =
% 
%      1     4     7
% 
% 
% ans =
% 
%      2     4     7
% 
% 
% ans =
% 
%      1     5     7
% 
% 
% ans =
% 
%      2     5     7
% 
% 
% ans =
% 
%      1     4     8
% 
% 
% ans =
% 
%      2     4     8
% 
% 
% ans =
% 
%      1     5     8
% 
% 
% ans =
% 
%      2     5     8
%
% Tom Rainforth 01/07/16

if numel(varargin)<=2
    s1 = numel(varargin{1});
    s2 = numel(varargin{2});
    E = cell(s1,s2);
    for n=1:s1
        for m=1:s2
            E{n,m} = [varargin{1}(n),varargin{2}(m)];
        end
    end
else
    new_in = varargin(1:end-1);
    E = enumerate_cells(new_in{:});
    sE = size(E);
    nEnd = numel(varargin{end});
    E = repmat(E,[ones(1,numel(sE)),nEnd]);
    nE = prod(sE);
    for n=1:(nE*nEnd)
        E{n} = [E{n},varargin{end}(ceil(n/nE))];
    end
end


