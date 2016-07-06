function h = moment_error_plotter(truth,var_name,block_length,b_invert,varargin)

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
                xlabel ''',''ylabel'',' '[num2str(block_length) '' poing moving average of squared error'']'];

block_size = [1,block_length];
            
for n=1:numel(varargin)
    sample_name = inputname(n+3);
    this_name_x = ['x_' sample_name];
    err_name = ['err_' sample_name];
    if b_invert
        eval([err_name ' = moving_average((truth(:,end:-1:1)-varargin{n}.' var_name '(:,end:-1:1)).^2, block_size);']);
    else
        eval([err_name ' = moving_average((truth-varargin{n}.' var_name ').^2, block_size);']);
    end
    plot_string = [plot_string ',' err_name];
end

plot_string = [plot_string, ');'];

eval(plot_string);