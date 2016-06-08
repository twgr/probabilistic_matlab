function C = empirical_covariance(samples,bIgnoreNaN,varargin)

n=1;
m=1;
while n<=numel(varargin)
    variables{m} = varargin{n}; %#ok<AGROW>
    if n==numel(varargin) || ischar(varargin{n+1})
        dimsToUse{m} = 1:size(samples.var.(variables{m}),2); %#ok<AGROW>
        n = n+1;
        m = m+1;
    else
        dimsToUse{m} = varargin{n}; %#ok<AGROW>
        n = n+2;
        m = m+1;
    end
end
        
if ~all(cellfun(@(x) isnumeric(samples.var.(x)), variables))
    error('Only valid for numeric variables');
end

if ~exist('bIgnoreNaN','var') || isempty(bIgnoreNaN)
    bIgnoreNaN = false;
elseif ischar(bIgnoreNaN)
    variables = [{bIgnoreNaN}, variables];
    bIgnoreNaN = false;
end

array_sizes_2 = cellfun(@numel, dimsToUse);
nX2 = sum(array_sizes_2);
ids_dim_2 = [0,cumsum(array_sizes_2)];

C = NaN(nX2,nX2);

if isfield(samples.var,'relative_weights')
    w = samples.var.relative_weights;
elseif isstruct(samples)
    w = ones(size(samples.var.(variables{1}),1),1);
else
    w = samples.relative_particle_weights;
end

for n_1 = 1:nX2
    for n_2 = n_1:nX2
        i_var_1 = find(n_1<=ids_dim_2,1)-1;
        ind_var_1 = n_1-(ids_dim_2(i_var_1));
        i_var_2 = find(n_2<=ids_dim_2,1)-1;
        ind_var_2 = n_2-(ids_dim_2(i_var_2));
                
        X1 = samples.var.(variables{i_var_1})(:,dimsToUse{n_1}(ind_var_1));
        X2 = samples.var.(variables{i_var_2})(:,dimsToUse{n_2}(ind_var_2));
        
        if issparse(X1)
            X1 = convert_to_full_array(X1);
            X2 = convert_to_full_array(X2);
        end
        
        w_local = w;
        
        if bIgnoreNaN
            bNaN = isnan(X1) | isnan(X2);
            w_local(bNaN) = 0;
            X1(bNaN) = 0;
            X2(bNaN) = 0;
        end
        
        mX1 = sum(X1.*w_local)/sum(w_local);
        mX2 = sum(X2.*w_local)/sum(w_local);
                
        sum_w = sum(w_local);
        scale = sum_w/(sum_w.^2-sum(w_local.^2));
        
        C(n_1,n_2) = scale*sum(w.*(X1-mX1).*(X2-mX2));
        C(n_2,n_1) = C(n_1,n_2);
    end
end

end
% 
% X = NaN(size(samples.var.(variables{1}),1),nX2);
% 
% ids_dim_2 = [0,cumsum(array_sizes_2)];
% 
% for n=1:numel(variables)
%     X(:,(ids_dim_2(n)+1):ids_dim_2(n+1)) = samples.var.(variables{n});
% end
% 
% if ~isfield(samples,'relative_weights') && isprop(samples,'relative_weights')
%     if bIgnoreNaN
%         bNaN = isnan(X);
%         X(bNaN) = 0;
%         m = mean(X,1);
%         X = bsxfun(@minus,X,m);
%         X(bNaN) = 0;
%         weights = double(~bNaN);
%     else
%         m = mean(X,1);
%         X = bsxfun(@minus,X,m);
%         weights = [];
%     end
% else
%     if bIgnoreNaN
%         bNaN = isnan(X);
%         X(bNaN) = 0;
%         m = sum(bsxfun(@times,samples.relative_weights,X),1)./...
%                     sum(bsxfun(@times,~bNaN,samples.relative_weights),1);
%         X = bsxfun(@minus,X,m);
%         X(bNaN) = 0;
%         weights = bsxfun(@times,samples.relative_weights,~bNaN);
%         m_w = sum(weights,1);
%         weights = size(X,1)*bsxfun(@rdivide,weights,m_w);
%         weights = sqrt(weights); % As these will appear as w^2 in the calcs.
%     else
%         m = 1/(sum(samples.relative_weights))*sum(...
%                     bsxfun(@times,samples.relative_weights,X),1);
%         X = bsxfun(@minus,X,m);
%         weights = repmat(sqrt(samples.relative_weights./mean(samples.relative_weights)),1,nX2);
%     end
% end
% 
% C = scaled_cov(X,weights);
%     
% end
% 
% function C = scaled_cov(Xc,weights)
% 
% if isempty(weights)
%     scale = size(Xc,1)-1;
% else
%     scale = weights'*weights-1;
%     scale(scale==0) = 1;
%     Xc = Xc.*weights;
% end
% C = (Xc' * Xc) ./scale; 
%   
% end