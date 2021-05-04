function res = signal(f,tri,signal_type,location_eval)
% signal(f,tri) computes the value of the signal given the value at 
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


 if (nargin == 2)
     location_eval = 'barycenter';
     signal_type = 'vertex';
 elseif (nargin == 3) && ~isempty(signal_type)
     location_eval = 'barycenter';
 elseif (nargin == 4) && isempty(signal_type)
     signal_type = 'vertex';
 end


if strcmpi(signal_type,'vertex')
    
    switch lower(location_eval)
        case 'barycenter'
            res = mean(reshape(f(tri),size(tri)),2);
            
        case 'middle' % middle of edges
            M = size(tri,2);
            if M==2
                res = (f(tri(:,1)) + f(tri(:,2)) ) /2;
            elseif M==3
                res = ( f(tri(:,[2,3,3])) + f(tri(:,[1,1,2]))  )/2;
            elseif M==1
		res = f(tri);
            end
            
        case 'vertices'
            res = reshape(f(tri),size(tri));
    end
    
elseif strcmpi(signal_type,'face')% if the signal given by user is constant on each face.
    % In that case length(f) == size(tri,1)

        switch lower(location_eval)
            case 'barycenter'
                 res = f;
        end

end





end
