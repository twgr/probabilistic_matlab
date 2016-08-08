function [X_1, other_output] = compose_two_sample_objects(X_1,X_2,i_assign,n_1,n_2,other_out_1,other_out_2,i_preorder_1,i_preorder_2)
%compose_two_sample_objects
%
% Combines two stack_object type objects according to an assignment order
% given by i_assign.  Workings are similar to struct_array_to_single_struct
%
% Required inputs:
%   X_1, X2 = Stack objects to combine.
%   i_assign = Array of reassignment ids.  Ids refer to output rather than
%              input
%   n_1, n_2 = Number of samples present in each object
%
% Optional inputs:
%   other_out_1, other_out_2 = Other variables to combine according to the
%              same ordering
%   i_preorder_1,i_preorder_2 = Allows reordering of the input variables.
%               This is necessary over / in addition to i_assign if
%               replicating variables.
%
% Outputs:
%   X_1 = First stack object after second has been added on
%   other_output = Combined other outputs
%
% Tom Rainforth 08/06/16

if ~exist('i_preorder_1','var')
    i_preorder_1 = [];
end

if ~exist('i_preorder_2','var')
    i_preorder_2= [];
end

fields_1 = fields(X_1.var);
fields_2 = fields(X_2.var);

variables = fields_1;
if ~all(strcmpi(fields_1,fields_2))  
    for m=1:numel(fields_2);
        if ~isempty(strcmpi(fields_1,fields_2{m}))
            variables = [variables; fields_2{m}]; %#ok<AGROW>
        end
    end
end

%variables = unique([fields_1; fields_2],'stable');

for n_v = 1:numel(variables)
    if any(strcmpi(fields_1,variables{n_v})) && ~isempty(i_preorder_1)
        X_1.var.(variables{n_v}) = X_1.var.(variables{n_v})(i_preorder_1,:);
    end
    
     if any(strcmpi(fields_2,variables{n_v})) && ~isempty(i_preorder_2)
        X_2.var.(variables{n_v}) = X_2.var.(variables{n_v})(i_preorder_2,:);
    end
    
    if n_1==0
        X_1.var.(variables{n_v})(i_assign,:) = X_1.var.(variables{n_v});
        continue
    elseif n_2==0
        X_1.var.(variables{n_v})(i_assign,:) = X_2.var.(variables{n_v});
        continue
    end
    
    if any(strcmpi(fields_1,variables{n_v})) && any(strcmpi(fields_2,variables{n_v}))
        if (size(X_1.var.(variables{n_v}),2)==size(X_2.var.(variables{n_v}),2)) && strcmp(class(X_1.var.(variables{n_v})),class(X_2.var.(variables{n_v})))
            if strcmpi(i_assign,':')
                X_1.var.(variables{n_v}) = [X_1.var.(variables{n_v});X_2.var.(variables{n_v})];
            else
                X_1.var.(variables{n_v})(i_assign,:) = [X_1.var.(variables{n_v});X_2.var.(variables{n_v})];
            end
        elseif isnumeric(X_1.var.(variables{n_v})) && isnumeric(X_2.var.(variables{n_v}))
            n_col_1 = size(X_1.var.(variables{n_v}),2);
            n_col_2 = size(X_2.var.(variables{n_v}),2);
            n_cols = max(n_col_1,n_col_2);
            X_1_mat = NaN(n_1,n_cols);
            X_1_mat(:,1:n_col_1) = X_1.var.(variables{n_v});
            X_2_mat = NaN(n_2,n_cols);
            X_2_mat(:,1:n_col_2) = X_2.var.(variables{n_v});
            if strcmpi(i_assign,':')
                X_1.var.(variables{n_v}) = [X_1_mat;X_2_mat];
            else
                X_1.var.(variables{n_v})(i_assign,:) = [X_1_mat;X_2_mat];
            end
        else
            % This should indicate that one is a cell and one is an array.
            if ~iscell(X_1.var.(variables{n_v}))
                X_1.var.(variables{n_v}) = num2cell(X_1.var.(variables{n_v}),2);
            end
            if ~iscell(X_2.var.(variables{n_v}))
                X_2.var.(variables{n_v}) = num2cell(X_2.var.(variables{n_v}),2);
            end
            if strcmpi(i_assign,':')
                X_1.var.(variables{n_v}) = [X_1.var.(variables{n_v});X_2.var.(variables{n_v})];
            else
                X_1.var.(variables{n_v})(i_assign,:) = [X_1.var.(variables{n_v});X_2.var.(variables{n_v})];
            end
        end
    elseif any(strcmpi(fields_1,variables{n_v}))
        if ~iscell(X_1.var.(variables{n_v}))
            X_1.var.(variables{n_v}) = num2cell(X_1.var.(variables{n_v}),2);
        end
        if strcmpi(i_assign,':')
            X_1.var.(variables{n_v}) = [X_1.var.(variables{n_v});cell(n_2,1)];
        else
            X_1.var.(variables{n_v})(i_assign,:) = [X_1.var.(variables{n_v});cell(n_2,1)];
        end
    elseif any(strcmpi(fields_2,variables{n_v}))
        if ~iscell(X_2.var.(variables{n_v}))
            X_2.var.(variables{n_v}) = num2cell(X_2.var.(variables{n_v}),2);
        end
        if strcmpi(i_assign,':')
            X_1.var.(variables{n_v}) = [cell(n_1,1);X_2.var.(variables{n_v})];
        else
            X_1.var.(variables{n_v})(i_assign,:) = [cell(n_1,1);X_2.var.(variables{n_v})];
        end
    else
        error('Should not be in variables in not in either true or false fields.  Somethings gone wrong in the code');
    end
end

if (~exist('other_out_1','var') || isempty(other_out_1)) && (~exist('other_out_2','var') || isempty(other_out_2))
    other_output = [];
else
    if (size(other_out_1,2)==size(other_out_2,2)) && strmpc(class(other_out_1),class(other_out_2))
        % These can be directly concatinated
    elseif (isnumeric(other_out_1) || isempty(other_out_1)) && (isnumeric(other_out_2) || isempty(other_out_2))
        n_col_1 = size(other_out_1,2);
        n_col_2 = size(other_out_2,2);
        n_cols = max(n_col_1,n_col_2);
        X_1_mat = NaN(n_1,n_cols);
        X_1_mat(1:size(other_out_1,1),1:n_col_1) = other_out_1;
        other_out_1 = X_1_mat;
        X_2_mat = NaN(n_2,n_cols);
        X_2_mat(1:size(other_out_2,1),1:n_col_2) = other_out_2;
        other_out_2 = X_2_mat;
    else
        % Cannot get out of num2cell in this scenario
        if ~iscell(other_out_1)
            other_out_1 = num2cell(other_out_1,2);
        end
        if ~iscell(other_out_2)
            other_out_2 = num2cell(other_out_2,2);
        end
        if isempty(other_out_1)
            other_out_1 = cell(n_1,size(other_out_2,2));
        end
        if isempty(other_out_2)
            other_out_2 = cell(n_2,size(other_out_1,2));
        end
        other_out_1 = [other_out_1, cell(size(other_out_1,1),size(other_out_2,2)-size(other_out_2,2))];
        other_out_2 = [other_out_2, cell(size(other_out_2,1),size(other_out_1,2)-size(other_out_2,2))];  
    end
    other_output(i_assign,:) = [other_out_1;other_out_2];
end

end