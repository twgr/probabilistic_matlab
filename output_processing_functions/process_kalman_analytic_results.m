A = dir;
A = A(3:end,:);
A = A(~[A(:).isdir],:);
names = {A(:).name};

error_names = {'Mean','Std_dev','Skewness','Excess_kurtosis'};

for n_name=1:numel(names)
    load(names{n_name});
    nB = regexp(names{n_name},'steps_')+6;
    if ~strcmpi(names{n_name}(nB+1),'.')
        nB = [nB,nB+1];
    end    
    for nx = 1:size(xsmooth,1)
        if size(xsmooth,1)==1
            name_var = 'x';
        else
            name_var = ['x_' num2str(nx)];
        end
        truths.(['b' names{n_name}(nB)]).(name_var).Mean = xsmooth(nx,:);
        truths.(['b' names{n_name}(nB)]).(name_var).Std_dev = sqrt(squeeze(Vsmooth(nx,nx,:))');
        truths.(['b' names{n_name}(nB)]).(name_var).Skewness = zeros(1,size(xsmooth,2));
        truths.(['b' names{n_name}(nB)]).(name_var).Excess_kurtosis = zeros(1,size(xsmooth,2));
    end
    truths.(['b' names{n_name}(nB)]).Cov_mat = Vsmooth;
    disp(n_name)
end