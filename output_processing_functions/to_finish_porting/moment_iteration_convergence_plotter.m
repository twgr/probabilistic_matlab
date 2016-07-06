function h = moment_iteration_convergence_plotter(truth,var_name,i_dims,varargin)

if ~isnumeric(i_dims)
    varargin = [{i_dims},varargin];
    i_dims = [];
    n_name = 2;
else
    n_name = 3;
end

if isempty(i_dims)
    i_dims = 1:size(truth,2);
end

plot_string = ['h = plot_distinct_lines(false,''title'',''' regexprep(var_name,'[_\.]',' ') ''',''xlabel'',''MCMC Iteration'',''ylabel'',''Sum of squared errors for statistic'',''Xlog'',true,''Ylog'',true'];
for n=1:numel(varargin)
    sample_name = inputname(n+n_name);
    this_name_x = ['x_' sample_name];
    err_name = ['err_' sample_name];
    eval([err_name ' = sum(bsxfun(@minus,truth(:,i_dims),varargin{n}.' var_name '(:,i_dims)).^2,2)/(numel(i_dims));']);
    plot_string = [plot_string ',' err_name];
end

plot_string = [plot_string, ');'];

eval(plot_string);