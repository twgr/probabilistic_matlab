function h = histogram_plotter(samples,field,ncols,nrows,d_vals,nbins,i_samples,ground_truths,xLims)

h = figure('units','normalized','outerposition',[0 0 1 1]);

if isfield(samples,'options')
    if isfield(samples.options,'iToKeep')
        d_vals = find(any(bsxfun(@eq,d_vals(:),samples.options.iToKeep(:)'),1));
    end
end

for n=1:nrows;
    for m=1:ncols;
        subplot(ncols,nrows,m+(n-1)*ncols);
        d = d_vals((n-1)*ncols+m);
        if d==0
            d = 1;
        end
        X = samples.var.(field)(:,d);
                
        if isfield(samples.var,'relative_weights') || isstruct(samples)
            % Support of outdated form
            w = samples.relative_weights;
        elseif isempty(samples.sparse_variable_relative_weights)
            if isempty(samples.relative_particle_weights)
                w = ones(size(X,1),1);
            else
                w = samples.relative_particle_weights;
            end
        elseif isnumeric(samples.sparse_variable_relative_weights)
            w = samples.sparse_variable_relative_weights(:,d);
        else
            w = samples.sparse_variable_relative_weights.(field)(:,d);
        end
        
        if exist('i_samples','var') && ~isempty(i_samples)
            
            
            n_particles = samples.options.n_particles;
            if ~strcmpi('smc',samples.options.inference_type) && ~strcmpi('msmc',samples.options.inference_type)
                n_iter = samples.options.n_iter;
                if isfield(samples.options,'n_nodes')
                    n_nodes = samples.options.n_nodes;
                else
                    n_nodes = 1;
                end
            else
                n_nodes = 1;
                if isfield(samples.options,'n_nodes')
                    n_iter = samples.options.n_nodes.*samples.options.n_iter;
                else
                    n_iter = samples.options.n_iter;
                end
                samples.options.n_nodes = n_nodes;
                samples.options.n_iter = n_iter;
            end
                        
            iteration_reorder = bsxfun(@plus,(1:n_particles)',(n_particles*n_iter)*(0:1:(n_nodes-1)));
            iteration_reorder = bsxfun(@plus,iteration_reorder(:),n_particles*(0:1:(n_iter-1)));
            iteration_reorder = iteration_reorder(:);
            
            X = X(iteration_reorder);
            w = w(iteration_reorder);
            
            if numel(i_samples)==1
                X = X(i_samples:end);
                w = w(i_samples:end);
            else
                X = X(i_samples);
                w = w(i_samples);
            end
        end
        
        if issparse(samples.var.(field))
            X = nonzeros(X);
            w = nonzeros(w);
        end
%         
%         if isfield(samples.var,'num_instances')
%             % Support of outdated form
%             [i,~,counts] = find(samples.var.num_instances(:,d));
%             w = counts.*samples.relative_weights(i);
%         elseif isfield(samples.var,'relative_weights') || isstruct(samples)
%             % Support of outdated form
%             w = samples.relative_weights;
%         elseif isempty(samples.sparse_variable_relative_weights)
%             if isempty(samples.relative_particle_weights)
%                 w = ones(size(X,1),1);
%             else
%                 w = samples.relative_particle_weights;
%             end
%         elseif isnumeric(samples.sparse_variable_relative_weights)
%             w = nonzeros(samples.sparse_variable_relative_weights(:,d));
%         else
%             w = nonzeros(samples.sparse_variable_relative_weights.(field)(:,d));
%         end
        
        w = w/sum(w);
            
        if ~exist('xLims','var') || isempty(xLims)        
            maxX = max(X);
            minX = min(X);
        else
            maxX = xLims((n-1)*ncols+m,2);
            minX = xLims((n-1)*ncols+m,1);
        end
        
        if maxX==minX
            plot([X(1);X(1)],[0;1],'b')
        else
            edges=(linspace(minX,maxX,nbins));
            [~, bin_assignment] = histc(full(X),edges);
            ibin_as = find(bin_assignment);
            counts = accumarray(bin_assignment(ibin_as),w(ibin_as),[numel(edges),1]);
            density = counts*nbins/(maxX-minX);
            
            %bar(0.5*(edges(2:end)+edges(1:end-1)),density,1)
            bar(edges,density,1);
        end
        
        if exist('ground_truths','var') && ~isempty(ground_truths)
            hold on;
            plot([ground_truths(d),ground_truths(d)],[0,max(density*1.2)],'r');
        end
        
        xlabel(['X_' num2str(d)]);
        
        if n==1 && m==1
            title(regexprep(inputname(1),'[_\.]',' '));
        end
        
        if n==1 && m==ncols
            title(['Probability densities of ' regexprep(field,'[_\.]',' ')]);
        end
    end
end

end