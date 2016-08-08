function squared_diffs = summary_squared_diffs(summary_1,summary_2)
%squared_diffs = summary_squared_diffs(summary_1,summary_2)
%
% Calculates squared differences between all fields of all variables for
% two summary objects (typically one that corresponds to results and one
% that is a ground truth).
%
% Inputs: 
%   summary_1, summary_2 = outputs from applying results_summary to a
%                          stack_object
% Outputs:
%   squared_diffs = Structure with a field with three fields.  Each has sub
%                   fields corresponding to common vields of summary_1 and
%                   summary_2.
%                       - all = All differences using bsxfun
%                       - pos_final = Differences at final iteration for
%                           each dimesion of the variable
%                       - conv_total = Difference summed along direction
%                           for each iteration.
%
% Tom Rainforth 08/08/16

var_names = fields(summary_1.var);

for n=1:numel(var_names)
    metric_names = fields(summary_1.var.(var_names{n}));
    for m=1:numel(metric_names)
        if isfield(summary_1.var.(var_names{n}),metric_names{m}) && ...
                isfield(summary_2.var.(var_names{n}),metric_names{m})
            squared_diffs.all.(var_names{n}).(metric_names{m}) = bsxfun(@minus,summary_1.var.(var_names{n}).(metric_names{m}),...
                                                                           summary_2.var.(var_names{n}).(metric_names{m})).^2;
            squared_diffs.pos_final.(var_names{n}).(metric_names{m}) = squared_diffs.all.(var_names{n}).(metric_names{m})(end,:);
            squared_diffs.conv_total.(var_names{n}).(metric_names{m}) = mean(squared_diffs.all.(var_names{n}).(metric_names{m}),2);
        end
    end
end

end