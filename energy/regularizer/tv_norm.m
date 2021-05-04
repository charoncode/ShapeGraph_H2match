% tv_norm.m: function that computes the total variation norm (TV norm) of a weight function (rho) defined on the edges 
% of a shape graph, and its gradient wrt the edge weights.
% Note: The gradient is taken wrt to a smooth approximation of the TV norm which uses the approximation |x| ≈ √(x + µ)^2 for small µ.
%
% Input:
%   shape: structure containing shape graph
%   tol  : small tolerance parameter for computing gradient (optional, default is 1e-8)
%
% Output:
%   tv   : TV norm of weight function rho
%   dtv  : gradient of the (approximate) TV norm with respect to the edge weights 

function [tv, dtv] = tv_norm(shape, tol)

% Set default tolerance to 1e-8
if nargin < 2
    tol = 1e-8;
end

% Extract weight function from shape graph structure
w = shape.rho; 
E = size(w,1); % number of edges

% Source and target edges of the edge adjacency matrix
[Is, It] = find(shape.E_adj); 
L_ind = find(Is < It);

% Compute TV norm
tv = sum(abs(w(It(L_ind)) - w(Is(L_ind))));

% Compute gradient of (approximate) TV norm
dtv = zeros(E,1);
for i = 1:length(L_ind)
   
   k = L_ind(i);
   
   % Stable version of TV flow
   dtv(Is(k)) = abs(w(Is(k)) - w(It(k))) * (w(Is(k)) - w(It(k)) ) / sqrt(tol^2 + (w(Is(k)) - w(It(k)) )^2 );
   dtv(It(k)) = abs(w(Is(k)) - w(It(k))) * (w(It(k)) - w(Is(k)) ) / sqrt(tol^2 + (w(It(k)) - w(Is(k)) )^2 );  
%  dtv(Is(k)) = (w(Is(k)) - w(It(k)) ) / sqrt(tol^2 + (w(Is(k)) - w(It(k)) )^2 );
%  dtv(It(k)) = (w(It(k)) - w(Is(k)) ) / sqrt(tol^2 + (w(It(k)) - w(Is(k)) )^2 );
   
end

end