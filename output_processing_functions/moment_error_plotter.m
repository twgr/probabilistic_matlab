function h = moment_error_plotter(truth,var_name,b_invert,varargin)

if isstruct(b_invert)
    varargin = [{b_invert}, varargin];
    b_invert = true;
end

if b_invert
    xlabel = 'Distance from end of chain';
else
    xlabel = 'Distance from start of chain';
end
    
plot_string = ['h = plot_distinct_lines(false,''title'',''' regexprep(var_name,'[_\.]',' ') ''',''xlabel'','''...
                xlabel ''',''ylabel'',''Sum of squared errors for statistic / ' xlabel ''''];
for n=1:numel(varargin)
    sample_name = inputname(n+2);
    this_name_x = ['x_' sample_name];
    err_name = ['err_' sample_name];
    if b_invert
        eval([err_name ' = cumsum((truth(:,end:-1:1)-varargin{n}.' var_name '(:,end:-1:1)).^2)./(1:size(truth,2));']);    
    else
        eval([err_name ' = cumsum((truth-varargin{n}.' var_name ')).^2./(1:size(truth,2));']); 
    end
    plot_string = [plot_string ',' err_name];
end

plot_string = [plot_string, ');'];

eval(plot_string);