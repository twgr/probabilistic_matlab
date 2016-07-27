function [x, iNonZeros] = get_variable_dim(samples,field,d,i_samples)
%get_variable_dim   Extracts a single variable dimension
%
% x = get_variable_dim(samples,field,d,i_samples)
%
% Variable extraction that removes all zeros when sparse.
%
% Inputs
%   samples = stack_object
%   field = Variable field
%   d = Dimension
%
% Outputs
%   x = Collapsed samples for that dimension
%   iNonZeros = Indices of non-zero points
%
% Tom Rainforth 27/07/16

x = samples.var.(field)(i_samples,d);
if issparse(x)
    x = nonzeros(x);
    if nargout>1
        iNonZeros = find(x(:,d));
    end
elseif nargout>1
    iNonZeros = 1:size(x_local,2);
end