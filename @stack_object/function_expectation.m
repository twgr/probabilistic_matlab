function mu = function_expectation(samples,f,b_supress_warning)
%varargout = function_expectation(samples,f,b_supress_warning)
%
% Calculates the expectation of a function f using the produced samples.
% 
% Inputs: samples,f
%   - f is the function to take an expectation of.  It should take in the
%     var struct from samples as input and return a column vector of
%     evaluations at each point or a matrix of multiple outputs where the
%     rows correspond to different samples
%
% Optional inputs: b_supress_warning
%   - Supresses the warning and need for human input in the recompression.
%     Use with caution as this can cause computer crashes if you cannot
%     store the uncompressed samples in memory!
%
% Outputs: mu = 1xNo row vector of means where No is the number of outputs
%               of f per sample
%
% Tom Rainforth 05/05/17

if ~isempty(samples.sparse_variable_relative_weights)
    warning('UNCOMP:WARN','Calculating function expectations requires uncompressing samples.\nConsider uncompressing upfront or providing f at the infer stage.');
    [~,msgId] = lastwarn;
    if ~isempty(msgId)
        warning('off',msgId);
    end
    samples = samples.uncompress_samples(b_supress_warning);
end

outs = f(samples.var);
w = samples.relative_particle_weights;
w = w/sum(w);
mu = sum(bsxfun(@times,outs,w),1);

end

