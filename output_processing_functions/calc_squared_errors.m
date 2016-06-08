function conv_summary = calc_squared_errors(conv_summary,truth_summary)

batch_names = fields(conv_summary);
var_names = fields(conv_summary.(batch_names{1}));
var_names = var_names(~cellfun(@(x) any(strcmpi(x,{'other_outputs','other_outpus','options','con'})), var_names));
error_names = fields(truth_summary.(batch_names{1}).(var_names{1}));

for n=1:numel(batch_names)
    for m=1:numel(var_names)
        for p=1:numel(error_names)
            conv_summary.(batch_names{n}).(var_names{m}).(error_names{p}) = bsxfun(@minus,...
                conv_summary.(batch_names{n}).(var_names{m}).(error_names{p}),...
                truth_summary.(batch_names{n}).(var_names{m}).(error_names{p})).^2;
            
            conv_summary.complete_state_sequence.(var_names{m}).(error_names{p})(:,n) = mean(conv_summary.(batch_names{n}).(var_names{m}).(error_names{p}),2);
            conv_summary.all_iter.(var_names{m}).(error_names{p})(n,:) = conv_summary.(batch_names{n}).(var_names{m}).(error_names{p})(end,:);
            conv_summary.(batch_names{n}).(var_names{m}) = rmfield(conv_summary.(batch_names{n}).(var_names{m}),(error_names{p}));
        end
    end
end

for n=1:numel(batch_names)
    for m=1:numel(var_names)
        conv_summary.(var_names{m}).ESS(:,n) = conv_summary.(batch_names{n}).(var_names{m}).ESS;
        conv_summary.(var_names{m}).ESS_first_iter(:,n) = conv_summary.(batch_names{n}).(var_names{m}).ESS_first_iter;
    end
end

end