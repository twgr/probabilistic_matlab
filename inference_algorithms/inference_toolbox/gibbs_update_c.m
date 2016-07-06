function [c, node_weights] = gibbs_update_c(log_Zs,c,n_cycles)
%gibbs_update_c Performs a gibbs update of conditional node indices
%
% Sweeps through the current conditional node indices and performs a Gibbs
% update whereby that conditional index can remain the same or switch to
% one nodes currently marked as unconditional.
%
% Required inputs:
%   log_Zs = Log marginal likelihood for each node as defined by eq 4 in
%            the paper
%   c      = Current conditional node indices
% Optional inputs:
%   n_cycles = Number of gibbs sweeps to perform, 1 by default.
%
% Outputs:
%   c      = Conditional node indices after update
%   node_weights = Normalized relative weights of the seperate nodes for
%                  Rao-Blackwellization
%
% Tom Rainforth 07/06/15

if ~exist('n_cycles','var') || isempty(n_cycles)
    % Number of gibbs sweep to perform.  Only 1 is required but the
    % accuracy of the node weights is improved by multiple sweeps
    n_cycles = 1;
end

P = numel(c);

% The scaling of the Zs does not matter as the node weights are only
% relative
Zs_unscaled = exp(log_Zs(:)-max(log_Zs));
zeta_sum = zeros(size(Zs_unscaled));

for nc = 1:n_cycles    
    for j = 1:P
        % Calculate zeta for this step
        zeta = Zs_unscaled;
        zeta(c(1:P~=j)) = 0;
        zeta = zeta/sum(zeta);       
        % Sample new c_j from categorical distribution
        cum_zeta = cumsum(zeta);
        cum_zeta(end) = 1;
        c_j = 1+sum(rand>cum_zeta);
        % Update in c vector
        c(j) = c_j;
        % Track the sum of all zeta for calculating the relative node
        % weights used in the Rao-Blackwellization of the samples
        zeta_sum = zeta_sum+zeta;
    end   
end

% The node weights are relative and so the absolute value is arbitrary.  We
% take the convention here that the weights sum to 1.  The weights as per
% the paper are P times the weights here
node_weights = zeta_sum/sum(zeta_sum);

end