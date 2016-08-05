function [x, iNonZeros] = get_variable(samples,field,dims,i_samples)
%get_variable   Extracts a variable
%
% x = get_variable(samples,field,dims,i_samples)
%
% Variable extraction that removes all zeros when sparse.
%
% Inputs
%   samples = stack_object
%   field = Variable field
%   dims = Dimensions
%
% Outputs
%   x = Collapsed samples for that dimension
%   iNonZeros = Indices of non-zero points
%
% Tom Rainforth 27/07/16

if ~exist('i_samples','var')
    i_samples = ':';
end

x = samples.var.(field)(i_samples,dims);
if issparse(x)
    x = nonzeros(x);
    if nargout>1
        iNonZeros = find(x);
    end
elseif nargout>1
    iNonZeros = (1:size(x,1))';
end