function Xfull = convert_to_full_array(Xsparse)
%convert_to_full_array
%
% Xfull = convert_to_full_array(Xsparse)
%
% Undoes the compression applied by the compress_samples function.  Note
% that this is different to just calling full because of encoding.
%
% Inputs:
%   Xsparse = Sparse array under compression encoding.  See
%             compress_samples for more details.
%
% Outputs:
%   Xfull = Full array uncoded
%
% Tom Rainforth 05/07/16

[iU,jU] = find(Xsparse);
Xfull = full(Xsparse);

for n=1:numel(iU)
    if n==numel(iU) || jU(n+1)~=jU(n)
        Xfull(iU(n)+1:end,jU(n)) = Xsparse(iU(n),jU(n));
    else
        Xfull((iU(n)+1):(iU(n+1)-1),jU(n)) = Xsparse(iU(n),jU(n));
    end
end
        