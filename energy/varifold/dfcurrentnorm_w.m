function [dxg,dfg,drhog]= dfcurrentnorm_w(fs1,fs2,objfun)
% DFCURRENTNORM(fs1,fs2,objfun) computes the derivative of the current 
% distance with respect to fs1.x and fs1.f.
%
% This function is equivalent to the function  "dfshape_kernel_distance" called with objfun.kernel_geom='gaussian', 
%    objfun.kernel_signal='gaussian' and objfun.kernel_grass='linear'. However, dfcurrentnorm
%    should be faster.
%
% Inputs:
%   fs1 : fshape structure containing the source
%   fs2 : fshape structure containing the target.
%   objfun is a struct with the parameters of the functional
%
% Outputs:
%   dxg is a (n xd) vector containing the gradient wrt fs1.x
%   dfg is a (n x1) vector containing the gradient wrt fs1.f
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2014)

[~,d]=size(fs1.x);
[Nx,~] = size(fs1.G);

sigmax=objfun.kernel_distance.kernel_size_geom;
sigmaf=objfun.kernel_distance.kernel_size_signal;

% Current corresponding to the fshapes templatef and target
[X,tf,Xix]=fcatoms(fs1.x,fs1.f,fs1.G,objfun.signal_type);
[Y,tg,Xiy]=fcatoms(fs2.x,fs2.f,fs2.G,objfun.data_signal_type);

WX=1+fs1.rho;
WY=1+fs2.rho;

%---------------------------------------------------------------------
% Computation of the gradients with respect to the X and F of the
% currents'representations
%--------------------------------------------------------------------

if (d<=2) && ~strcmp(objfun.kernel_distance.method,'matlab') % zeros padding
    X = [X,zeros(size(X,1),1)];
    Xix = [Xix,zeros(size(Xix,1),1)];
    Y = [Y,zeros(size(Y,1),1)];
    Xiy = [Xiy,zeros(size(Xiy,1),1)];
end

% switch objfun.kernel_distance.method
%     case 'cuda' %Mex/Cuda files with exact computations of exponential
%         
%         DXig=2*(GaussFunGpuConv([X,tf]',[X,tf]',Xix',sigmax,sigmaf)-GaussFunGpuConv([X,tf]',[Y,tg]',Xiy',sigmax,sigmaf))';
%         
%         AAA=GaussFunGpuConvgrad([X,tf]',[X,tf]',Xix',sigmax,sigmaf);
%         BBB=GaussFunGpuConvgrad([X,tf]',[Y,tg]',Xiy',sigmax,sigmaf);
%         DXg=2*[sum(Xix' .* (AAA(:,:,1)-BBB(:,:,1)),1);sum(Xix' .* (AAA(:,:,2)-BBB(:,:,2)),1);sum(Xix' .* (AAA(:,:,3) -BBB(:,:,3)),1)]'; % valid only for d=3
%         
%         Dtfg=2*sum(Xix .* (GaussFunGpuConvgradfun([X tf]',[X tf]',Xix',sigmax,sigmaf)-GaussFunGpuConvgradfun([X tf]',[Y tg]',Xiy',sigmax,sigmaf))',2);
%         
%     case 'mexc'
% 
%             DXig=2*(dsOptim([X,tf],[X,tf],[sigmax sigmaf],Xix,0)-dsOptim([X tf],[Y tg],[sigmax sigmaf],Xiy,0));
%             
%             AAA=dsOptim([X,tf],[X,tf],[sigmax sigmaf],Xix,1);
%             BBB=dsOptim([X,tf],[Y,tg],[sigmax sigmaf],Xiy,1);
%             DXg=2*[sum(Xix .* (AAA(:,:,1)-BBB(:,:,1)),2),sum(Xix .* (AAA(:,:,2)-BBB(:,:,2)),2),sum(Xix .* (AAA(:,:,3) -BBB(:,:,3)),2)];
%             
%             Dtfg=2*sum(Xix .* (dsOptim([X tf],[X tf],[sigmax sigmaf],Xix,2)-dsOptim([X tf],[Y tg],[sigmax sigmaf],Xiy,2)),2);
%         
%     otherwise  %pure Matlab version
        Ny=size(Y,1);
        
        % Calcul de A=exp(-|x_i -x_j|^2/(2*lam^2))
        
        Sx=zeros(Nx);Sxy=zeros(Nx,Ny);
        Dxx=zeros(Nx);Dxy=zeros(Nx,Ny);
        for l=1:d
            Sx=Sx+(repmat(X(:,l),1,Nx)-repmat(X(:,l)',Nx,1)).^2;
            Sxy=Sxy+(repmat(X(:,l),1,Ny)-repmat(Y(:,l)',Nx,1)).^2;

            Dxx=Dxx+(repmat(Xix(:,l),1,Nx).*repmat(Xix(:,l)',Nx,1));
            Dxy=Dxy+(repmat(Xix(:,l),1,Ny).*repmat(Xiy(:,l)',Nx,1));
        end
        
        Aff=rho((repmat(tf,1,Nx)-repmat(tf',Nx,1)).^2,0,sigmaf);
        Afg=rho((repmat(tf,1,Ny)-repmat(tg',Nx,1)).^2,0,sigmaf);
        Apff=rho((repmat(tf,1,Nx)-repmat(tf',Nx,1)).^2,1,sigmaf);
        Apfg=rho((repmat(tf,1,Ny)-repmat(tg',Nx,1)).^2,1,sigmaf);
        
        Wxx=WX*WX';
        Wxy=WX*WY';  
        
        Axx=rho(Sx,0,sigmax);
        Apxx=rho(Sx,1,sigmax).*Dxx.*Aff.*Wxx;
        Axy=rho(Sxy,0,sigmax);
        Apxy=rho(Sxy,1,sigmax).*Dxy.*Afg.*Wxy;
              
        
        DXg=zeros(Nx,d);
        
        for l=1:d
            DXg(:,l)=4*( sum((Apxx.*(repmat(X(:,l),1,Nx)-repmat(X(:,l)',Nx,1))),2)-sum(Apxy.*(repmat(X(:,l),1,Ny)-repmat(Y(:,l)',Nx,1)),2) ); % scalar kernel case
        end
        
        DXig=2*((Axx.*Aff.*Wxx)*Xix-(Axy.*Afg.*Wxy)*Xiy); % dot product part
        Dtfg=4*(sum(Axx.*Dxx.*Apff.*Wxx.*(repmat(tf,1,Nx)-repmat(tf',Nx,1)),2)-sum(Axy.*Dxy.*Apfg.*Wxy.*(repmat(tf,1,Ny)-repmat(tg',Nx,1)),2));
        
        drhog = 2*(sum(Axx.*Dxx.*Aff.*repmat(WX,1,Nx),2)-sum(Axy.*Dxy.*Afg.*repmat(WY',Nx,1),2));
%end


[dxg,dfg]= chain_rule(fs1.x,fs1.G,DXg,DXig,Dtfg,objfun.signal_type);

end

