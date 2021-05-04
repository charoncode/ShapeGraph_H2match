function res = dsignal(df,tri,nb_points,type_signal)
% DSIGNAL(df,tri) computes the value of the derivative of the signal. given the value at 
% the vertices (ie barycentric coordinate)
%
% Input:
%  f: signal value at each vertices (n x 1 column vector)
%  tri: list of edges (T x M matrix of indices)
%  type_signal : optional string. If  type_signal=='face', f is assumed to be
%  the values of the signal at the center of the faces else f is assumed to
%  be the values of the signal at the vertices (ie pts)
%
% output:
%  f: signal value at the center of the cells (n x 1 column vector)
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, G. Nardi, A. Trouve (2012-2014)

if nargin==3 || isempty(type_signal)
    type_signal = 'barycenter';
end
[~,M] = size(tri);

switch lower(type_signal)
    case 'face'
        % if the signal given by user is constant on each face. 
            res = df(:,1);
    otherwise
        % if the signal given by user is the value at each vertex. 
        res = accumarray(tri(:),repmat(df(:,1),M,1),[nb_points,1],[],0)/M;


end




end