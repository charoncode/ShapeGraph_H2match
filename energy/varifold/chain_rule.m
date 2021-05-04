function [dxg2,dfg2]=chain_rule(x,tri,DXg,DXig,Dtfg,type_signal)
%
% Computation of the gradient with respect to x by distributing the previous gradients
% on points of the initial shape. (Chain's rule). There is two versions of the same
% code somewhat equivalent : one using matlab for loop (that is not so bad
% due to JIT precompilation process) and one using built in matlab function
% accumarray.
%
% Input:
%  x: list of vertices (n x d matrix)
%  G: list of edges (T x M matrix)
%  DXg: derivative wrt X (center of the cells: T x d matrix)
%  DXig: derivative wrt Xi (p-vectors: T x d matrix)
%  Dtfg: derivative wrt tf (signal at center of the cells: T x1 colunm vector)
%
% Output:
%  dxg: derivative wrt x (n x d matrix)
%  dfg: derivative wrt f (n x 1 colunm vector)
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2014)


[nx,d]=size(x);
[~,M] = size(tri);

if M==1 % point cloud case  chain's rule stops here
    dfg2 = Dtfg;
    dxg2 = DXg;
    return
end

dfg2 = dsignal(Dtfg,tri,nx,type_signal);

if (d==2)
    U2 =  [accumarray(tri(:),repmat(DXg(:,1),M,1),[nx,1],[],0),...
           accumarray(tri(:),repmat(DXg(:,2),M,1),[nx,1],[],0)] / M;
elseif (d==3)
    U2 =  [accumarray(tri(:),repmat(DXg(:,1),M,1),[nx,1],[],0),...
           accumarray(tri(:),repmat(DXg(:,2),M,1),[nx,1],[],0),...
           accumarray(tri(:),repmat(DXg(:,3),M,1),[nx,1],[],0)] / M;
end


if (M==2) % curve case
    if (d==2)
        dxg2 = U2 + [accumarray(tri(:),[-DXig(:,1);DXig(:,1)] ,[nx,1],[],0),...
                     accumarray(tri(:),[-DXig(:,2);DXig(:,2)] ,[nx,1],[],0) ];
    elseif (d==3)
        dxg2 = U2 + [accumarray(tri(:),[-DXig(:,1);DXig(:,1)] ,[nx,1],[],0),...
                     accumarray(tri(:),[-DXig(:,2);DXig(:,2)] ,[nx,1],[],0),...
                     accumarray(tri(:),[-DXig(:,3);DXig(:,3)] ,[nx,1],[],0) ];
    end
    
elseif (M==3) && (d==3) % surface case
    
    Xa=x(tri(:,1),:);
    Xb=x(tri(:,2),:);
    Xc=x(tri(:,3),:);
    

    [dU1,dU2,dU3,dV1,dV2,dV3] = dcross((Xb-Xa),(Xc-Xa),DXig/2);
    
    dxg2 = U2 + [accumarray(tri(:),[-dU1-dV1;+dU1 ;+dV1 ],[nx,1],[],0),...
                 accumarray(tri(:),[-dU2-dV2;+dU2 ;+dV2 ],[nx,1],[],0),...
                 accumarray(tri(:),[-dU3-dV3;+dU3 ;+dV3 ],[nx,1],[],0)];

elseif (M==4) && (d==3) % simplexes case

    Xa=x(tri(:,1),:);
    Xb=x(tri(:,2),:);
    Xc=x(tri(:,3),:);
    Xd=x(tri(:,4),:);
    
    u=(Xb-Xa);v=(Xc-Xa);w=(Xd-Xa);

    %[dU1,dU2,dU3,dV1,dV2,dV3,dW1,dW2,dW3] =  dcrosss(u,v,w,H);
    H = DXig/6;
    
    dU1 = ( v(:,2) .* w(:,3) - v(:,3) .* w(:,2) ).* H;
    dU2 = ( v(:,3) .* w(:,1) - v(:,1) .* w(:,3) ).* H;
    dU3 = ( v(:,1) .* w(:,2) - v(:,2) .* w(:,1) ).* H;
   
    dV1 = ( w(:,2) .* u(:,3) - w(:,3) .* u(:,2) ).* H;
    dV2 = ( w(:,3) .* u(:,1) - w(:,1) .* u(:,3) ).* H;
    dV3 = ( w(:,1) .* u(:,2) - w(:,2) .* u(:,1) ).* H;
    
    dW1 = ( u(:,2) .* v(:,3) - u(:,3) .* v(:,2) ).* H;
    dW2 = ( u(:,3) .* v(:,1) - u(:,1) .* v(:,3) ).* H;
    dW3 = ( u(:,1) .* v(:,2) - u(:,2) .* v(:,1) ).* H;
    

    
    dxg2 = U2 + [accumarray(tri(:),[-dU1-dV1-dW1 ;+dU1 ;+dV1 ; dW1],[nx,1],[],0),...
                  accumarray(tri(:),[-dU2-dV2-dW2 ;+dU2 ;+dV2 ; dW2],[nx,1],[],0),...
                  accumarray(tri(:),[-dU3-dV3-dW3 ;+dU3 ;+dV3 ; dW3],[nx,1],[],0)];
    

    
elseif (M==3) && (d==2) % simplexes case

    Xa=x(tri(:,1),:);
    Xb=x(tri(:,2),:);
    Xc=x(tri(:,3),:);
    
    U=(Xb-Xa);V=(Xc-Xa);
    
    H = DXig/2;
    
    dU1 =   V(:,2) .* H;
    dU2 =  -V(:,1) .* H;
    
    dV1 =  -U(:,2) .* H ;
    dV2 =   U(:,1) .* H;
    
    dxg2 = U2 + [accumarray(tri(:),[-dU1-dV1;+dU1 ;+dV1 ],[nx,1],[],0)...
                 accumarray(tri(:),[-dU2-dV2;+dU2 ;+dV2 ],[nx,1],[],0)];
    
end

end


function [dxg,dfg]= Chain_rule1(x,G,DXg,DXig,Dtfg) %#ok<DEFNU>
% Equivalent code with loop... but almost no speedup due to JIT compilation process in last matlab versions.
[nx,d]=size(x);
[~,M] = size(G);

if M==1 % (measure case), chain's rule stops here
    dfg = Dtfg;
    dxg = DXg;
    return
end

% if M==2 (curves case)|| M==3 (surface)
% the case of simplexes is not implemented yet (ie M == d+1)
U=zeros(nx,d);
dfg=zeros(nx,1);

for l=1:d
    for m=1:M
        U(:,l)=U(:,l)+nmaps(DXg(:,l),G(:,m),nx);
    end
end
U=U/M; % because X_i=(x_i1+..+x_iM)/M

for m=1:M
    dfg=dfg+nmaps(Dtfg,G(:,m),nx);
end
dfg=dfg/M; % because F_i=(f_i1+..+f_iM)/M

if M==2 % curve case
    
    for l=1:d
        U(:,l) = U(:,l)+nmaps(DXig(:,l)  ,G(:,2),nx) ;
        U(:,l) = U(:,l)-nmaps(DXig(:,l)  ,G(:,1),nx);
    end
    dxg=U;
    
elseif M==3 && d==3 % surface case
    
    Xa=x(G(:,1),:);
    Xb=x(G(:,2),:);
    Xc=x(G(:,3),:);
    
    u=Xb-Xa; v=Xc-Xa;
    
    % Automatic differentiation (still manual however...)
    dU=v(:,3).*DXig(:,1)-v(:,1).*DXig(:,3);
    U(:,2)=U(:,2)+nmaps(dU,G(:,2),nx)-nmaps(dU,G(:,1),nx);
    
    dU=-v(:,2).*DXig(:,1)+v(:,1).*DXig(:,2);
    U(:,3)=U(:,3)+nmaps(dU,G(:,2),nx)-nmaps(dU,G(:,1),nx);
    
    
    dU=-v(:,3).*DXig(:,2)+v(:,2).*DXig(:,3);
    U(:,1)=U(:,1)+nmaps(dU,G(:,2),nx)-nmaps(dU,G(:,1),nx);
    
    
    dU=-u(:,3).*DXig(:,1)+u(:,1).*DXig(:,3);
    U(:,2)=U(:,2)+nmaps(dU,G(:,3),nx)-nmaps(dU,G(:,1),nx);
    
    dU=+u(:,2).*DXig(:,1)-u(:,1).*DXig(:,2);
    U(:,3)=U(:,3)+nmaps(dU,G(:,3),nx)-nmaps(dU,G(:,1),nx);
    
    dU=+u(:,3).*DXig(:,2)-u(:,2).*DXig(:,3);
    U(:,1)=U(:,1)+nmaps(dU,G(:,3),nx)-nmaps(dU,G(:,1),nx);
    
    dxg=U;
    
end
end

function V = nmaps(U,C,m)
% maps
%   This piece of code maps the value of U on the vecor V of size m through
%   the allocation vector C ie V(i) goes at location C(i) in U. C is not
%   suppose to be one to one and several values of V maps at the same
%   location are just added. More precisely, this code implements the
%   following code using spare matrices for speed optimization.
%
% Equivalent of the built-in Matlab function "accumarray"... but for loop are efficient with JIT pre compilation.
%
% See also : accumarray

V=zeros(m,1);
for i=1:length(C)
    V(C(i))=V(C(i))+U(i);
end

end
