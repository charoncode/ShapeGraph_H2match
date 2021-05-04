% huber.m: function to computes the Huber function and its derivative.
%
% Input:
%  v            : point at which we evaluate the Huber function 
%                (Note: since the huber function is a function from R to R, v should in 
%                 theory be a scalar. But v can be a vector, in which case we evaluate
%                 the Huber function element-wise at each entry of v)
%  lambda, gamma: constants for the huber function

% Output:
%   hub  : value of huber function evaluated at v
%   dhub : derivative of the huber function evaluated at v

function [hub, dhub] = huber(v, lambda, gamma)

% Dimensions
M = length(v);

% Pre allocate arrays
hub = zeros(M,1);
dhub = zeros(M,1);

% Define constant
q = lambda/gamma;

% Find entries where |v| <= q and where |v| > q
ind_less = find(abs(v) <= q); 
ind_more = find(abs(v) > q);

% Compute Huber function and derivative for each case
hub(ind_less) = 0.5*gamma*v(ind_less).^2;
dhub(ind_less) = gamma*v(ind_less);

hub(ind_more) = lambda*(abs(v(ind_more)) - q/2);
dhub(ind_more) = lambda*sign(v(ind_more));

end


