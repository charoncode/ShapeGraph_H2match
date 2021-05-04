function res = pVectors(pts,tri,scheme)
% computes (a representation of) the p-vector (ie tangent vector for curve and normals for surfaces). 
%
% Input:
%  pts: list of vertices (n x d matrix)
%  tri: list of edges (T x M matrix of indices)
%  scheme: discretization scheme, 'forward' or 'centered' (for curves only)
%
% Output:
%  res: list of p-vectors (d x M matrix)
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2014)

if nargin==2
    scheme='forward';
end

M = size(tri,2);
d = size(pts,2);
N = size(pts,1);
res = zeros(N,d);

if (M==2) % curves
    
    
    if strcmp(scheme,'forward')
        
        res=(pts(tri(:,2),:)-pts(tri(:,1),:));
        
    else
        for i=1:d
            res(:,i) = gradient(pts(:,i),1);
        end
            
%         if 0 %norm(pts(end,:)-pts(1,:))<10^(-6) || isequal(tri(end,:),[N 1])
%             
%             res(1,:)=(pts(2,:)-pts(end,:))/(2*N);
%             res(end,:)=(pts(1,:)-pts(end-1,:))/(2*N);
%             
%         end
               
    end
    
    
elseif (M==3) && (d==3) % surfaces

    u=pts(tri(:,2),:)-pts(tri(:,1),:); 
    v=pts(tri(:,3),:)-pts(tri(:,1),:);
    
    res =[u(:,2).*v(:,3)-u(:,3).*v(:,2),...
          u(:,3).*v(:,1)-u(:,1).*v(:,3),...
          u(:,1).*v(:,2)-u(:,2).*v(:,1)];

elseif (M==1) % points

    res = repmat( 1./size(tri,1) ,size(tri,1),1)  ;
  
elseif (M==3) && (d==2)% simplexes

    u=pts(tri(:,2),:)-pts(tri(:,1),:); 
    v=pts(tri(:,3),:)-pts(tri(:,1),:);

    res =  u(:,1).*v(:,2) - u(:,2).*v(:,1);

elseif (M==4) && (d==3)% simplexes

    u=pts(tri(:,2),:)-pts(tri(:,1),:); 
    v=pts(tri(:,3),:)-pts(tri(:,1),:);
    w=pts(tri(:,4),:)-pts(tri(:,1),:);

    res =  u(:,1).*v(:,2).*w(:,3) + v(:,1).*w(:,2).*u(:,3) + w(:,1).*u(:,2).*v(:,3)...
         - u(:,3).*v(:,2).*w(:,1) - v(:,3).*w(:,2).*u(:,1) - w(:,3).*u(:,2).*v(:,1); %det with Sarrus formula

end

end
