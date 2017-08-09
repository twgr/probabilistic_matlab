function [i_resample,n_times_sampled] = resample_indices(w,n_samples,method)
%resample_indices   Samples the indices for resampling
%
% Given weights, a number of samples and a method, calculates the indices
% of the resampled particles.
%
% Inputs:
%   w (vector of +ve reals) = Weights, need not be normalized
%   n_samples (+ve integer) = Number of samples to take
%   method ('stratified' | 'systematic' | 'multinomial' | 'residual')
%       = Employed resampling method, Default = 'systematic'
%
% Outputs
%   i_resample = Indices of the particles to keep
%   n_times_sampled = Number of times each particle was sampled
%
% Tom Rainforth 13/06/16

if ~exist('method','var') || isempty(method)
    method = 'systematic';
end

% Normalize weights
w = w/sum(w);
if strcmpi(method,'residual')
    Nw = n_samples*w;
    base_n = floor(Nw);
    R = sum(base_n);
    residual = Nw-base_n;
    i_base = [];
    base_n_dec = base_n;
    while (numel(i_base)<R)
        i_base = [i_base;find(base_n_dec)]; %#ok<AGROW>
        base_n_dec = max(0,base_n_dec-1);
    end
    if R==n_samples
        n_times_sampled = base_n;
        i_resample = i_base;
    else
        % Stratified resampling on the rest
        [i_res,n_res] = resample_indices(residual,n_samples-R,'stratified');
        i_resample = [i_base;i_res];
        n_times_sampled = base_n+n_res;
    end
else
    switch method
        case 'multinomial'
            drawsForResample = rand(n_samples,1);
        case 'stratified'
            drawsForResample = rand(n_samples,1)/n_samples+(0:(n_samples-1))'/n_samples;
        case 'systematic'
            drawsForResample = rand/n_samples+(0:(n_samples-1))'/n_samples;
        otherwise
            error('Unrecognized resampling method');
    end
    edges = min([0;cumsum(w)],1);
    edges(end) = 1;
    [n_times_sampled,i_resample] = histc(drawsForResample,edges);
    n_times_sampled = n_times_sampled(1:end-1);
    i_resample = i_resample(:);
    n_times_sampled = n_times_sampled(:);
end