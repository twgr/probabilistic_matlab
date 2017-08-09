function samples = uncompress_samples(samples,b_supress_warning)

if ~exist('b_supress_warning','var') || isempty(b_supress_warning)
    b_supress_warning = false;
end

S = whos('samples');
if S.bytes>1e7 && ~b_supress_warning
    % Realistic chance of swamping memory
    button = questdlg(sprintf('Uncompressing can crash your computer if not enough memory\npresent. Curernt size of %d MB may increase significantly\nContinue?',...
                               S.bytes/1e6),'WARNING!','Yes','No','No');
    if ~strcmpi(button,'Yes')
        return;
    end
end

if numel(samples)~=1
    for n=1:numel(samples)
        samples(n) = uncompress_samples(samples(n),true);
    end
    return
end

p_fields = fields(samples.var);
for n_f=1:numel(p_fields)
    samples.var.(p_fields{n_f}) = convert_to_full_array(samples.var.(p_fields{n_f}));
end

samples.sparse_variable_relative_weights = [];