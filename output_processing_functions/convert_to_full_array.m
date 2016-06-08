function Xfull = convert_to_full_array(Xsparse)

[iU,jU] = find(Xsparse);
Xfull = full(Xsparse);

for n=1:numel(iU)
    if n==numel(iU) || jU(n+1)~=jU(n)
        Xfull(iU(n)+1:end,jU(n)) = Xsparse(iU(n),jU(n));
    else
        Xfull((iU(n)+1):(iU(n+1)-1),jU(n)) = Xsparse(iU(n),jU(n));
    end
end
        