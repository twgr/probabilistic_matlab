function [h,hsubs,densities,edges] = histogram_plotter(samples,field,n_bins,...
            nrows_plot,ncols_plot,b_discrete,dims,i_samples,edges,ground_truths)
%histogram_plotter    Used to plot generate a setup histogram subplots
%
% [h,densities,edges] = histogram_plotter(samples,field,n_bins,...
%     nrows_plot,ncols_plot,b_discrete,dims,i_samples,edges,ground_truths)
%
% Call results_binning and the plots the results as subplot on a single
% figure.
% Inputs:
%   samples,filed,n_bins = See results_binning.m
%   n_rows_plot,n_cols_plot = Defines the subplot array dimensions
%   b_discrete,dims,i_samples,edges = See results_binning.m
%   ground_truths = Vector of "true" values that will added to each plot
%                   as vertical lines
% Outputs:
%   h = Handle for overall figure object
%   hsubs = Vector of handles for the subplots
%   densities = See results_binning.m
%   edges = See results_binning.m
%
% Tom Rainforth 27/07/16
        
[~,n_dims_total] = size(samples.var.(field));

if ~exist('nrows_plot','var') || isempty(nrows_plot)
    nrows_plot = 1;
end

if ~exist('ncols_plot','var') || isempty(ncols_plot)
    ncols_plot = 1;
end

if ~exist('b_discrete','var')
    b_discrete = [];
end

if ~exist('dims','var') || isempty(dims)
    dims = 1:n_dims_total;
end

if ~exist('i_samples','var')
    i_samples = [];
end

if ~exist('edges','var')
    edges = [];
end

[densities,edges,b_distinct_values] = results_binning(samples,field,n_bins,b_discrete,dims,i_samples,edges);

h = figure('units','normalized','outerposition',[0 0 1 1]);
hsubs = NaN(size(dims));

for n=1:nrows_plot;
    for m=1:ncols_plot;
        nd = (n-1)*ncols_plot+m;
        try
            d = dims(nd);
        catch
            continue
        end
                
        subplot(ncols_plot,nrows_plot,nd);
        if ~b_distinct_values(nd)
            hsubs(nd) = plot([edges(1,nd);edges(1,nd)],[0;1],'b');
        else
            hsubs(nd) = bar(0.5*(edges(1:end-1,nd)+edges(2:end,nd)),densities(:,nd),1);
        end
        
        if exist('ground_truths','var') && ~isempty(ground_truths)
            hold on;
            plot([ground_truths(nd),ground_truths(nd)],[0,max(densities(:,nd))*1.2],'r');
        end
        
        xlabel(['X_{' num2str(d) '}']);
%         
%         if n==1 && m==1
%             title(regexprep(inputname(1),'[_\.]',' '));
%         end
%         
        if n==1 && m==1
            title(['Probability densities of ' regexprep(field,'[_\.]',' ') ' for ' regexprep(inputname(1),'[_\.]',' ')]);
        end
        
        set(gca,'FontSize',12);
        
    end
end

end