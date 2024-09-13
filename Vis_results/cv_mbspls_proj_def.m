function Xs = cv_mbspls_proj_def(Xs, weights)
%
%   Perform projection deflation of matrices X and Y using u and v
%   Please check Monteiro et al. 2016 for details:
%   doi:10.1016/j.jneumeth.2016.06.011
%
%   Inputs: X, Y    - data matrices in the form: samples x features.
%
%           u, v    - weight vectors for X and Y, respectively
%
%
%   Outputs: X, Y    - deflated data matrices
%
%
%   Version: 2016-07-06
%__________________________________________________________________________

% Written by Clara Vetter

for num_m=1:size(Xs,2)
    Xs{num_m} = Xs{num_m} - (Xs{num_m}*weights{num_m})*weights{num_m}';
end