function [dxg,dfg,drhog]= dfvarifoldnorm_binet_w(fs1,fs2,objfun)
% FVARIFOLDNORM_BINET(fs1,fs2,objfun) computes  derivative of the fvarifold 
% distance with Cauchy-Binet kernel with respect to fs1.x and
% fs1.f
%
% This function is equivalent to the function  "dfshape_kernel_distance" called with objfun.kernel_geom='gaussian', 
%    objfun.kernel_signal='gaussian' and objfun.kernel_grass='binet'. However, dfvarifold_binet
%    should be faster.
%
% Inputs:
%   fs1: fshape structure containing the source
%   fs2: fshape structure containing the target.
%   objfun: is a structure containing the data attachment term parameters
%
% Outputs:
%   dxg is a (nd x1) vector containing the gradient of g wrt x
%   dfg is a (n x1) vector containing the gradient of g wrt f
%
% See also : dfvarifold_gauss, dfshape_kernel_distance, dmatchterm
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2014)

[~,d]=size(fs1.x);
Nx = size(fs1.G,1);

sigmax=objfun.kernel_distance.kernel_size_geom;
sigmaf=objfun.kernel_distance.kernel_size_signal;

% Current corresponding to the fshapes templatef and target
[X,tf,Xix]=fcatoms(fs1.x,fs1.f,fs1.G,objfun.signal_type);
[Y,tg,Xiy]=fcatoms(fs2.x,fs2.f,fs2.G,objfun.data_signal_type);

WX=1+fs1.rho;
WY=1+fs2.rho;

% Filter out degenerate faces
M_areaX=sqrt(sum(Xix.^2,2));
ind=find(M_areaX<1e-5*mean(M_areaX));
% X=X(ind,:);
% tf=tf(ind,:);
% Xix=Xix(ind,:);

%---------------------------------------------------------------------
% Computation of the gradients with respect to the X and F of the
% currents'representations
%--------------------------------------------------------------------

if (d<=2) && ~strcmp(objfun.kernel_distance.method,'matlab')  % zeros padding
    d=3; % can be messy...
    X = [X,zeros(size(X,1),1)];
    Xix = [Xix,zeros(size(Xix,1),1)];
    Y = [Y,zeros(size(Y,1),1)];
    Xiy = [Xiy,zeros(size(Xiy,1),1)];
end

sqrtnormXix = sqrt(sqrt(sum(Xix .^2,2)));
Xixtilde = Xix ./  repmat(sqrtnormXix,1,size(Xix,2));
sqrtnormXiy = sqrt(sqrt(sum(Xiy .^2,2)));
Xiytilde = Xiy ./  repmat(sqrtnormXiy,1,size(Xiy,2));

Ny=size(Y,1);

% switch objfun.kernel_distance.method
%     case 'cuda' %Mex/Cuda files with exact computations of exponential
%         DXitildeg = zeros(Nx,d);
%         for l=1:d
%             DXitildegxx = sum(GaussFunGpuConv([X tf]',[X tf]',(Xixtilde.* repmat( Xixtilde(:,l),1,3))',sigmax,sigmaf)' .* Xixtilde,2);
%             DXitildegxy = sum(GaussFunGpuConv([X tf]',[Y tg]',(Xiytilde.* repmat( Xiytilde(:,l),1,3))',sigmax,sigmaf)' .* Xixtilde,2);
%             
%             DXitildeg(:,l) = 4 * ( DXitildegxx- DXitildegxy);
%         end
%         
%         % Expend the squared scalar product
%         Xiixtilde = Xixtilde.^2;
%         Xiiytilde = Xiytilde.^2;
%         Xijxtilde = sqrt(2) * Xixtilde .* Xixtilde(:,[3 1 2]);
%         Xijytilde = sqrt(2) * Xiytilde .* Xiytilde(:,[3 1 2]);
%         
%         
%         
%         XX = GaussFunGpuConvgradfun([X tf]',[X tf]',Xiixtilde',sigmax,sigmaf)';
%         Xx = GaussFunGpuConvgradfun([X tf]',[X tf]',Xijxtilde',sigmax,sigmaf)';
%         XY =  GaussFunGpuConvgradfun([X tf]',[Y tg]',Xiiytilde',sigmax,sigmaf)';
%         Xy =  GaussFunGpuConvgradfun([X tf]',[Y tg]',Xijytilde',sigmax,sigmaf)';
%         
%         Dtfg = 2 * ( sum(Xiixtilde .* (XX - XY),2) +  sum(Xijxtilde .*(Xx - Xy) ,2));
%         
%         XX = GaussFunGpuConvgrad([X tf]',[X tf]',Xiixtilde',sigmax,sigmaf);
%         Xx = GaussFunGpuConvgrad([X tf]',[X tf]',Xijxtilde',sigmax,sigmaf);
%         XY =  GaussFunGpuConvgrad([X tf]',[Y tg]',Xiiytilde',sigmax,sigmaf);
%         Xy =  GaussFunGpuConvgrad([X tf]',[Y tg]',Xijytilde',sigmax,sigmaf);
%         
%         DXg =2* [sum(Xiixtilde .* (XX(:,:,1)-XY(:,:,1))'+  Xijxtilde .* (Xx(:,:,1)-Xy(:,:,1))',2),...
%             sum(Xiixtilde .* (XX(:,:,2)-XY(:,:,2))'+  Xijxtilde .* (Xx(:,:,2)-Xy(:,:,2))',2),...
%             sum(Xiixtilde .* (XX(:,:,3)-XY(:,:,3))'+  Xijxtilde .* (Xx(:,:,3)-Xy(:,:,3))',2)];
%         
%             
%     otherwise
        
        % Calcul de A=exp(-|x_i -x_j|^2/(2*lam^2))
        Sx=zeros(Nx);
        Sxy=zeros(Nx,Ny);
        Dxx=zeros(Nx);
        Dxy=zeros(Nx,Ny);
        
        for l=1:d
            Sx=Sx+(repmat(X(:,l),1,Nx)-repmat(X(:,l)',Nx,1)).^2;
            Sxy=Sxy+(repmat(X(:,l),1,Ny)-repmat(Y(:,l)',Nx,1)).^2;

            Dxx=Dxx+(repmat(Xixtilde(:,l),1,Nx).*repmat(Xixtilde(:,l)',Nx,1));
            Dxy=Dxy+(repmat(Xixtilde(:,l),1,Ny).*repmat(Xiytilde(:,l)',Nx,1));
        end
        
        Aff=rho((repmat(tf,1,Nx)-repmat(tf',Nx,1)).^2,0,sigmaf);
        Afg=rho((repmat(tf,1,Ny)-repmat(tg',Nx,1)).^2,0,sigmaf);
        Apff=rho((repmat(tf,1,Nx)-repmat(tf',Nx,1)).^2,1,sigmaf);
        Apfg=rho((repmat(tf,1,Ny)-repmat(tg',Nx,1)).^2,1,sigmaf);
        
        Axx=rho(Sx,0,sigmax); %AxxFF=Axx.*Aff;
        Apxx=rho(Sx,1,sigmax);%.*Dxx.^2 .*Aff;
        Axy=rho(Sxy,0,sigmax); %AxyFG=Axy.*Afg;
        Apxy=rho(Sxy,1,sigmax);%.*Dxy.^2 .*Afg;
        
        Wxx=WX*WX';
        Wxy=WX*WY';
        
        DXg=zeros(Nx,d);
        
        for l=1:d
            DXg(:,l)=4*(sum((Apxx.*Dxx.^2 .*Aff.*Wxx.*(repmat(X(:,l),1,Nx)-repmat(X(:,l)',Nx,1))),2)...
                -sum(Apxy.*Dxy.^2 .*Afg.*Wxy.*(repmat(X(:,l),1,Ny)-repmat(Y(:,l)',Nx,1)),2)); % scalar kernel case
        end
        
        DXitildeg=4*((Axx.*Aff.*Dxx)*Xixtilde-(Axy.*Afg.*Dxy)*Xiytilde); % dot product part
        
        Dtfg=4*(sum(Axx.*Dxx.^2 .*Apff.*Wxx.*(repmat(tf,1,Nx)-repmat(tf',Nx,1)),2)-sum(Axy.*Dxy.^2 .*Apfg.*Wxy.*(repmat(tf,1,Ny)-repmat(tg',Nx,1)),2));
        
        drhog = 2*(sum(Axx.*Dxx.^2.*Aff.*repmat(WX,1,Nx),2)-sum(Axy.*Dxy.^2.*Afg.*repmat(WY',Nx,1),2));
        
%end

% Computation of the gradient with respect to x by distributing the previous gradients on points of the initial shape.
% Filter out degenerate faces
DXg(ind,:)=0;
DXitildeg(ind,:)=0;

DXig=zeros(Nx,d);
    
if d==2
    for i = 1:Nx
        % Jacobian matrix D_\xi \tilde\xi
        DXixitilde = [Xix(i,1).^2 + 2*Xix(i,2).^2,-Xix(i,1).* Xix(i,2);...
                     -Xix(i,1).* Xix(i,2)        , Xix(i,2).^2+ 2*Xix(i,1).^2 ]./ (2 * sqrtnormXix(i) .^5);
        DXig(i,:) =  (DXixitilde * DXitildeg(i,:)')';
    end
elseif d==3
    for i = 1:Nx
        % Jacobian matrix D_\xi \tilde\xi
        DXixitilde = [Xix(i,1).^2 + 2*Xix(i,2).^2 + 2*Xix(i,3).^2,-Xix(i,1).* Xix(i,2),-Xix(i,1).* Xix(i,3);...
                     -Xix(i,2).* Xix(i,1)       , Xix(i,2).^2+ 2*Xix(i,1).^2+ 2*Xix(i,3).^2,-Xix(i,2).* Xix(i,3)  ;...
                     -Xix(i,3).* Xix(i,1)       ,-Xix(i,3).* Xix(i,2) , Xix(i,3).^2+ 2*Xix(i,1).^2+ 2*Xix(i,2).^2]./ (2 * sqrtnormXix(i) .^5);
        DXig(i,:) =  (DXixitilde * DXitildeg(i,:)')';
    end 
end

% End of the chain's rule
[dxg,dfg]= chain_rule(fs1.x,fs1.G,DXg,DXig,Dtfg,objfun.signal_type);

end
