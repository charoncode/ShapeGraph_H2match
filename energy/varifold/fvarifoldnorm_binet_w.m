function g= fvarifoldnorm_binet_w(fs1,fs2,objfun)
% FVARIFOLDNORM_BINET(fs1,fs2,objfun) computes the functional varifold
% distance with Cauchy-Binet kernel. Possible method are
% 'cuda', 'mexc' or 'matlab'.
%
% This function is equivalent to the function  "fshape_kernel_distance" called with objfun.kernel_geom='gaussian', 
%    objfun.kernel_signal='gaussian' and objfun.kernel_grass='binet'. However, fvarifold_binet
%    should be faster.
%
% Inputs:
%   fs1 : fshape structure containing the source
%   fs2 : fshape structure containing the target.
%   objfun : is a structure containing the data attachment term parameters
%
% Output
%   g : a real number.
%
% See also : fvarifold_gauss, fshape_kernel_distance, matchterm
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2014)

d=size(fs1.x,2);
M=size(fs1.G,2);

sigmax=objfun.kernel_distance.kernel_size_geom;
sigmaf=objfun.kernel_distance.kernel_size_signal;

[X,tf,Xix]=fcatoms(fs1.x,fs1.f,fs1.G,objfun.signal_type);
[Y,tg,Xiy]=fcatoms(fs2.x,fs2.f,fs2.G,objfun.data_signal_type);

WX=1+fs1.rho;
WY=1+fs2.rho;

sqrtnormXix = sqrt(sqrt(sum(Xix .^2,2)));
Xix = Xix ./  repmat(sqrtnormXix,1,size(Xix,2));
sqrtnormXiy = sqrt(sqrt(sum(Xiy .^2,2)));
Xiy = Xiy ./  repmat(sqrtnormXiy,1,size(Xiy,2));


if (d<=2) && ~strcmp(objfun.kernel_distance.method,'matlab') % zeros padding
    X = [X,zeros(size(X,1),1)];
    Xix = [Xix,zeros(size(Xix,1),1)];
    Y = [Y,zeros(size(Y,1),1)];
    Xiy = [Xiy,zeros(size(Xiy,1),1)];
end

% switch objfun.kernel_distance.method
%     case 'cuda'  
%         
%         Xiix = Xix.^2;
%         Xiiy = Xiy.^2;
%         Xijx = sqrt(2) * Xix .* Xix(:,[3 1 2]); % only true if d==3
%         Xijy = sqrt(2) * Xiy .* Xiy(:,[3 1 2]); % only true if d==3
%         
%         % norm(x)=
%         XX= GaussFunGpuConv([X,tf]',[X,tf]',Xiix',sigmax,sigmaf)';
%         Xx= GaussFunGpuConv([X,tf]',[X,tf]',Xijx',sigmax,sigmaf)';
%         Nx = Xiix(:)' * XX(:) +  Xijx(:)' * Xx(:);
%                
%         % morm(y) =
%         YY= GaussFunGpuConv([Y,tg]',[Y,tg]',Xiiy',sigmax,sigmaf)';
%         Yy= GaussFunGpuConv([Y,tg]',[Y,tg]',Xijy',sigmax,sigmaf)';
%         Ny = Xiiy(:)' * YY(:) + Xijy(:)' * Yy(:) ;
%         
%         %prs(x,y) =
%         XY= GaussFunGpuConv([X,tf]',[Y,tg]',Xiiy',sigmax,sigmaf)';
%         Xy= GaussFunGpuConv([X,tf]',[Y,tg]',Xijy',sigmax,sigmaf)';
%         Pxy =  Xiix(:)' * XY(:) + Xijx(:)' *Xy(:);
%         
%         g= Nx + Ny - 2 * Pxy;
%     case 'mexc' % use mex file dsOptim
%             
%             Xiix = Xix.^2;
%             Xiiy = Xiy.^2;
%             Xijx = sqrt(2) * Xix .* Xix(:,[3 1 2]); % only true if d==3
%             Xijy = sqrt(2) * Xiy .* Xiy(:,[3 1 2]); % only true if d==3
%             
%             % norm(x)=
%             XX= dsOptim([X,tf],[X,tf],[sigmax,sigmaf],Xiix,0);
%             Xx= dsOptim([X,tf],[X,tf],[sigmax,sigmaf],Xijx,0);
%             Nx = Xiix(:)' * XX(:) +  Xijx(:)' * Xx(:);
%             
%             
%             % morm(y) =
%             YY= dsOptim([Y,tg],[Y,tg],[sigmax,sigmaf],Xiiy,0);
%             Yy= dsOptim([Y,tg],[Y,tg],[sigmax,sigmaf],Xijy,0);
%             Ny = Xiiy(:)' * YY(:) + Xijy(:)' * Yy(:) ;
%             
%             %prs(x,y) =
%             XY= dsOptim([X,tf],[Y,tg],[sigmax,sigmaf],Xiiy,0);
%             Xy= dsOptim([X,tf],[Y,tg],[sigmax,sigmaf],Xijy,0);
%             Pxy =  Xiix(:)' * XY(:) + Xijx(:)' *Xy(:);
%             
%             g= Nx + Ny - 2 * Pxy;
% 
%     otherwise  %pure matlab version
        
        Nx=size(X,1);
        Ny=size(Y,1);
        
        % Calcul de A=exp(-|x_i -x_j|^2/(2*lam^2))
        Sx=zeros(Nx);
        Dxx=zeros(Nx);
        Sy=zeros(Ny);
        Dyy=zeros(Ny);
        Sxy=zeros(Nx,Ny);
        Dxy=zeros(Nx,Ny);
        
        
        for l=1:d
            Sx=Sx+(repmat(X(:,l),1,Nx)-repmat(X(:,l)',Nx,1)).^2;
            Sy=Sy+(repmat(Y(:,l),1,Ny)-repmat(Y(:,l)',Ny,1)).^2;
            Sxy=Sxy+(repmat(X(:,l),1,Ny)-repmat(Y(:,l)',Nx,1)).^2;

            Dxx=Dxx+ (repmat(Xix(:,l),1,Nx).*repmat(Xix(:,l)',Nx,1));
            Dyy=Dyy+(repmat(Xiy(:,l),1,Ny).*repmat(Xiy(:,l)',Ny,1));
            Dxy=Dxy+(repmat(Xix(:,l),1,Ny).*repmat(Xiy(:,l)',Nx,1));
        end
        
        Ax=rho(Sx,0,sigmax);
        Ay=rho(Sy,0,sigmax);
        Axy=rho(Sxy,0,sigmax);
        
        %par.kernel_size_mom=par.objfun.kernel_size_signal;
        FF=rho((repmat(tf,1,Nx)-repmat(tf',Nx,1)).^2,0,sigmaf);
        GG=rho((repmat(tg,1,Ny)-repmat(tg',Ny,1)).^2,0,sigmaf);
        FG=rho((repmat(tf,1,Ny)-repmat(tg',Nx,1)).^2,0,sigmaf);
        
        WXX=WX*WX';
        WYY=WY*WY';
        WXY=WX*WY';
       
        g=sum(sum(Ax.*Dxx.^2 .*FF.*WXX))+sum(sum(Ay.*Dyy.^2 .*GG.*WYY))-2*sum(sum(Axy.*Dxy.^2 .*FG.*WXY));
%end

end



