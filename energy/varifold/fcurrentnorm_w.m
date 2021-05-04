function g= fcurrentnorm_w(fs1,fs2,objfun)
% FCURRENTNORM(fs1,fs2,objfun) computes (the square) od the radial functional 
% currents norm between two fshapes fs1 and fs2.
%
% This function is equivalent to the function  "fshape_kernel_distance" called with objfun.kernel_geom='gaussian', 
%    objfun.kernel_signal='gaussian' and objfun.kernel_grass='linear'. However, fcurrentnorm
%    should be faster.
%
% Inputs:
%   fs1 : fshape structure containing the source
%   fs2 : fshape structure containing the target.
%   objfun : is a structure containing the data attachment term parameters
% Output
%   g : a real number.
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2014)

d=size(fs1.x,2);
M=size(fs1.G,2);

sigmax=objfun.kernel_distance.kernel_size_geom;
sigmaf=objfun.kernel_distance.kernel_size_signal;

[X,tf,Xix]=fcatoms_w(fs1.x,fs1.f,fs1.G,fs1.rho,objfun.signal_type);
[Y,tg,Xiy]=fcatoms_w(fs2.x,fs2.f,fs2.G,fs2.rho,objfun.data_signal_type);

Xix=Xix.*repmat(1+fs1.rho,1,size(Xix,2));
Xiy=Xiy.*repmat(1+fs2.rho,1,size(Xiy,2));

if M<=2 && (strcmp('mexc',objfun.kernel_distance.method))
    objfun.kernel_distance.method = 'matlab';
    warning('Switch to matlab version in matchterm as M<=2.')
end


% switch objfun.kernel_distance.method
%     case 'cuda' % use mex/cuda file  with brut force double summing
%         
%         % norm(x)=
%         XX = GaussFunGpuConv([X,tf]',[X,tf]',Xix',sigmax,sigmaf)';
%         Nx = Xix(:)' * XX(:);
%         
%         
%         % morm(y) =
%         YY = GaussFunGpuConv([Y,tg]',[Y,tg]',Xiy',sigmax,sigmaf)';
%         Ny = Xiy(:)' * YY(:);
%         
%         
%         %prs(x,y) =
%         XY= GaussFunGpuConv([X,tf]',[Y,tg]',Xiy',sigmax,sigmaf)';
%         Pxy =  Xix(:)' * XY(:);
%         
%         
%         g= Nx + Ny - 2 *  Pxy;
%         
%     case 'mexc'
%         
%             % norm(x) =
%             XX = dsOptim([X,tf],[X,tf],[sigmax sigmaf],Xix,0);
%             Nx = Xix(:)' * XX(:);
%             
%             % morm(y) =
%             YY= dsOptim([Y,tg],[Y,tg],[sigmax sigmaf],Xiy,0);
%             Ny = Xiy(:)' * YY(:);
%             
%             %prs(x,y) =
%             XY= dsOptim([X,tf],[Y,tg],[sigmax sigmaf],Xiy,0);
%             Pxy =  Xix(:)' * XY(:);
%             
%             g= Nx + Ny - 2 *  Pxy ;
% 
%     otherwise
        Nx=size(X,1);
        Ny=size(Y,1);
        
        % Calcul de A=exp(-|x_i -x_j|^2/(2*lam^2))
        Sx=zeros(Nx);
        Sy=zeros(Ny);
        Sxy=zeros(Nx,Ny);        
        Dxx=zeros(Nx);
        Dyy=zeros(Ny);
        Dxy=zeros(Nx,Ny);
        
        for l=1:d
            Sx=Sx+(repmat(X(:,l),1,Nx)-repmat(X(:,l)',Nx,1)).^2;
            Sy=Sy+(repmat(Y(:,l),1,Ny)-repmat(Y(:,l)',Ny,1)).^2;
            Sxy=Sxy+(repmat(X(:,l),1,Ny)-repmat(Y(:,l)',Nx,1)).^2;

            Dyy=Dyy+(repmat(Xiy(:,l),1,Ny).*repmat(Xiy(:,l)',Ny,1));          
            Dxx=Dxx+(repmat(Xix(:,l),1,Nx).*repmat(Xix(:,l)',Nx,1));
            Dxy=Dxy+(repmat(Xix(:,l),1,Ny).*repmat(Xiy(:,l)',Nx,1));
        end
        
        
        
        Ax=rho(Sx,0,sigmax);
        Ay=rho(Sy,0,sigmax);
        Axy=rho(Sxy,0,sigmax);
        
        %par.kernel_size_mom=par.objfun.kernel_size_signal;
        FF=rho((repmat(tf,1,Nx)-repmat(tf',Nx,1)).^2,0,sigmaf);
        GG=rho((repmat(tg,1,Ny)-repmat(tg',Ny,1)).^2,0,sigmaf);
        FG=rho((repmat(tf,1,Ny)-repmat(tg',Nx,1)).^2,0,sigmaf);
        
        % Norm(x)^2 + Norm(y)^2 - 2 * prs(X,Y)
        g=sum(sum(Ax.*Dxx.*FF))+sum(sum(Ay.*Dyy.*GG))-2*sum(sum(Axy.*Dxy.*FG));

%end

end