function [s, other_outputs] = struct_array_to_single_struct(s_array,assignment_order,other_outputs_array)
%struct_array_to_single_struct
%
% Converts an array of structures (or objects) to a single structure /
% object according to an assignment_order.  This is complicated by the
% fact that variables may change type or dimension between the different
% instances in the array due to the dynanic typing.  If both are numeric,
% padding with NaNs is used, otherwise different instances are stored as a
% cell array in output.
%
% Inputs:
%  s_array = Starting structure / object array
%  assignment_order = Cell array dictating the assignment order.  These
%                      are indices on the outputs such that
%                  s.field(assignment_order{1}(1)) = s_array(1).field(1)
%  other_outputs_array = Cell array of same size to be rearranged to a
%                        single array (either cell or numeric depending on
%                        type) using the same assignment_order
%
% Output:
%  s = Single structure / object with desired rearrangement
%  other_outputs = Rearrange other outputs in appropriate form
%
% Tom Rainforth 08/06/16

if ~exist('assignment_order','var')
    assignment_order = [];
end

name_that_wont_be_used = 'asfdh4myrq34qfwqr234hansd';

if exist('other_outputs_array','var')
    for n=1:numel(other_outputs_array)
        s_array(n).(name_that_wont_be_used) = other_outputs_array{n};
    end
end

variables = fields(s_array);

for n_v = 1:numel(variables)
    this_var = {s_array(:).(variables{n_v})}';
    if all(cellfun(@isempty,this_var))
        continue
    end
    
    sizes = cellfun(@(x) size(x,2), this_var);
    classes = cellfun(@class, this_var, 'UniformOutput', false);
    
    size_and_class = cellfun(@(x,y) [x,y],arrayfun(@num2str,sizes,'UniformOutput', false),classes,'UniformOutput',false);
    [unique_sizes_and_classes,~,type_belonging] = unique(size_and_class);
    n_unique = size(unique_sizes_and_classes,1);
    
    if n_unique==1
        if ~isempty(assignment_order)
            s.(variables{n_v})(cell2mat(assignment_order),:) = cell2mat(this_var);
        else
             s.(variables{n_v}) = cell2mat(this_var);
        end
    else
        temp_stack_objects = repmat(stack_object,n_unique,1);   
        if isempty(assignment_order)
            io = (1:numel(type_belonging))';
        else
            io = assignment_order;
        end
        for n_u = 1:n_unique
            temp_stack_objects(n_u).var.temp = cell2mat(this_var(type_belonging==n_u,:));
            temp_stack_objects(n_u).var.ids = cell2mat(io(type_belonging==n_u));
        end     
        iter_check = 1;
        while n_unique>1 && iter_check<100           
            if mod(n_unique,2)==1
                i_stack = randi(n_unique-1);
                temp_stack_objects(i_stack) = compose_temps(temp_stack_objects(i_stack),temp_stack_objects(end));
                temp_stack_objects = temp_stack_objects(1:end-1);
            end
            
            for n_s = 1:round(n_unique/2)
               temp_stack_objects(n_s) = compose_temps(temp_stack_objects(n_s),temp_stack_objects(end));
               temp_stack_objects = temp_stack_objects(1:end-1);
            end
            n_unique = round(n_unique/2);
            iter_check = iter_check+1;
            disp(iter_check);
        end
        s.(variables{n_v})(temp_stack_objects.var.ids,:) = temp_stack_objects.var.temp;
    end
end

if isfield(s,name_that_wont_be_used)
    other_outputs = s.(name_that_wont_be_used);
    s = rmfield(s,name_that_wont_be_used);
else
    other_outputs = [];
end

end

function temp = compose_temps(temp_obj_1,temp_obj_2)
    
   X_temp = stack_object;
   i_1 = (1:size(temp_obj_1.var.temp,1))';
   i_2 = size(temp_obj_1.var.temp,1)+(1:size(temp_obj_2.var.temp,1))';
   i_assign = [i_1;i_2];
   n_1 = size(i_1,1);
   n_2 = size(i_2,1);
   
   temp = compose_two_sample_objects(X_temp,temp_obj_1,temp_obj_2,i_assign,n_1,n_2,[],[]);
end
