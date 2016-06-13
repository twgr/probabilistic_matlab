function [samples,b_compress] = memory_check(samples,n_to_take,n_sampling_functions)
%memory_check  Checks if going to cause memory issues
%
% Used to turn on compression automatically if appear to be in danger of
% causing memory issues.
%
% Inputs
%   samples = Current samples generated
%   n_to_take = The number of times the provided samples will be repeated
%               upon completion
%   n_sampling_functions = Total number of sampling functions
%
% Outputs
%   samples = Same as input but compressed if required
%   b_compress = Whether compression should now be turned on
%
% Tom Rainforth 13/06/16

S = whos('samples');
s_mem = S.bytes*n_to_take;
if s_mem>5e7
    try
        memory_stats = memory;
        largest_array = memory_stats.MaxPossibleArrayBytes;
    catch
        % memory function is only availible in windows
        largest_array = 4e9;
    end
    
    if S.bytes*n_total > (largest_array/20)
        warning('In danger of swamping memory and crashing, turning b_compress on');
        b_compress = true;
        samples = compress_samples(samples, n_sampling_functions);
    end
end