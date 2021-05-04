function g= fshape_kernel_distance_w(fs1,fs2,objfun)
% FSHAPE_KERNEL_DISTANCE(templatefinal,target,objfun) computes kernel based
% distances between fshapes.
%
%  \sum_i\sum_j K_signal(f_i-g_j)^2) K_geom(-norm(x_i-y_j)^2) K_tan(angle(V_i,W_j))
% 
% Possible method are 'cuda' or 'matlab'.
%
% Inputs:
%   templatefinal : "fshape structure" containing the shooted template
%   target : "fshape structure" containing the target.
%   objfun : is a structure containing the data attachment term parameters (mandatory fields are : 'kernel_geom', 'kernel_signal' and 'kernel_grass' and the associated bandwidth)
% Output
%   g : a real number.
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2016)

d=size(fs1.x,2);
m=size(fs1.G,2)-1;

% discretize the fshapes
[center_faceX,signalX,normalsX]=fcatoms(fs1.x,fs1.f,fs1.G,objfun.signal_type);
[center_faceY,signalY,normalsY]=fcatoms(fs2.x,fs2.f,fs2.G,objfun.data_signal_type);

WX=1+fs1.rho;
WY=1+fs2.rho;

options = objfun.kernel_distance;

if min(m, d-m) ==0 % for points clouds or simplexes dim or codim == 0 : some simplifications occurs 

     g = ptcloud_or_simplexes_kernel_distance(center_faceX,signalX,normalsX,WX,center_faceY,signalY,normalsY,WY,options);

elseif min(m,d-m) ==1 % for curves or surface dim or codim ==1;

     % add a correcting weight if the target has missing cells
     if isfield(options,'weight_missing_data')
        normalsY = normalsY / options.weight_missing_data;
     end
     g = curves_or_surfaces_kernel_distance(center_faceX,signalX,normalsX,WX,center_faceY,signalY,normalsY,WY,options);

end

end



function g = ptcloud_or_simplexes_kernel_distance(center_faceX,signalX,normalsX,WX,center_faceY,signalY,normalsY,WY,objfun)
% this function is equivalent to the curves_or_surfaces_kernel_distance
% below with a radial_function_sphere constant to 1. This version avoid
% some computation and will be faster for fshape of dim or codim == 0.
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2016)


%switch objfun.method
%      case 'cuda'  % use cuda to speedup the computation
% 
% 	 eval(['fshape_scp=@fsimplex_scp_Gpu_',lower(objfun.kernel_geom),lower(objfun.kernel_signal),';']);
% 
%          % norm(x)^2 =
% 	 XX= fshape_scp(center_faceX',center_faceX',signalX',signalX',normalsX',normalsX',objfun.kernel_size_geom,objfun.kernel_size_signal);
% 	 PXX =  sum(XX);
%                   
%          % morm(y)^2 =
%          YY= fshape_scp(center_faceY',center_faceY',signalY',signalY',normalsY',normalsY',objfun.kernel_size_geom,objfun.kernel_size_signal);
% 	 PYY =  sum(YY);
%          
%          %prs(x,y) =
% 	 XY= fshape_scp(center_faceX',center_faceY',signalX',signalY',normalsX',normalsY',objfun.kernel_size_geom,objfun.kernel_size_signal);
% 	 PXY =  sum(XY);

%     otherwise

        [Nx,d]=size(center_faceX);
        Ny=size(center_faceY,1);
        
        %compute squared distances and angles
        distance_signalXX = (repmat(signalX,1,Nx)-repmat(signalX',Nx,1)).^2;
        distance_signalYY = (repmat(signalY,1,Ny)-repmat(signalY',Ny,1)).^2;
        distance_signalXY = (repmat(signalX,1,Ny)-repmat(signalY',Nx,1)).^2;
        
        distance_center_faceXX = zeros(Nx);
        distance_center_faceYY=zeros(Ny);
        distance_center_faceXY=zeros(Nx,Ny);
    
        for l=1:d
            distance_center_faceXX = distance_center_faceXX+(repmat(center_faceX(:,l),1,Nx)-repmat(center_faceX(:,l)',Nx,1)).^2;
            distance_center_faceYY = distance_center_faceYY+(repmat(center_faceY(:,l),1,Ny)-repmat(center_faceY(:,l)',Ny,1)).^2;
            distance_center_faceXY = distance_center_faceXY+(repmat(center_faceX(:,l),1,Ny)-repmat(center_faceY(:,l)',Nx,1)).^2;
        end
        
        % Geometric kernel      
        Kernel_geomXX = radial_function_geom(distance_center_faceXX,0,objfun);
        Kernel_geomYY = radial_function_geom(distance_center_faceYY,0,objfun);
        Kernel_geomXY = radial_function_geom(distance_center_faceXY,0,objfun);
        
        % Signal kernel
        Kernel_signalXX = radial_function_signal(distance_signalXX,0,objfun);
        Kernel_signalYY = radial_function_signal(distance_signalYY,0,objfun);
        Kernel_signalXY = radial_function_signal(distance_signalXY,0,objfun);

        % Area
        totalWX=normalsX.*repmat(WX,1,size(normalsX,2));
        totalWY=normalsY.*repmat(WY,1,size(normalsY,2));
        AreaXX = (totalWX * totalWX');
        AreaYY = (totalWY * totalWY');
        AreaXY = (totalWX * totalWY');
        
        % norm(x)=
        PXX = sum(sum(AreaXX .* Kernel_geomXX .* Kernel_signalXX));
        
        % morm(y) =
        PYY =sum(sum(AreaYY .* Kernel_geomYY .* Kernel_signalYY));
        
        %prs(x,y) =
        PXY = sum(sum(AreaXY .* Kernel_geomXY .* Kernel_signalXY));
        
%end         

 g= PXX + PYY - 2* PXY;
end






function g = curves_or_surfaces_kernel_distance(center_faceX,signalX,normalsX,WX,center_faceY,signalY,normalsY,WY,objfun)
%
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2016)

%switch objfun.method
%      case 'cuda'  % use cuda to speedup the computation
% 
% 	 eval(['fshape_scp=@fshape_scp_Gpu_',lower(objfun.kernel_geom),lower(objfun.kernel_signal),lower(objfun.kernel_grass),';']);
% 
%          % norm(x)^2 =
%          XX= fshape_scp(center_faceX',center_faceX',signalX',signalX',normalsX',normalsX',objfun.kernel_size_geom,objfun.kernel_size_signal,objfun.kernel_size_grass);
%          PXX =  sum(XX);
%          
%          
%          % morm(y)^2 =
%          YY= fshape_scp(center_faceY',center_faceY',signalY',signalY',normalsY',normalsY',objfun.kernel_size_geom,objfun.kernel_size_signal,objfun.kernel_size_grass);
%          PYY =  sum(YY);
%          
%          %prs(x,y) =
%          XY= fshape_scp(center_faceX',center_faceY',signalX',signalY',normalsX',normalsY',objfun.kernel_size_geom,objfun.kernel_size_signal,objfun.kernel_size_grass);
%          PXY =  sum(XY);
%                   
%      otherwise

        [Nx,d]=size(center_faceX);
        Ny=size(center_faceY,1);

        % Compute unit normals
        norm_normalsX = sqrt(sum(normalsX .^2,2));
        norm_normalsY = sqrt(sum(normalsY .^2,2));
        
        unit_normalsX = normalsX ./  repmat(norm_normalsX,1,size(normalsX,2));
        unit_normalsY = normalsY ./  repmat(norm_normalsY,1,size(normalsY,2));
        
        %compute squared distances and angles
        distance_signalXX = (repmat(signalX,1,Nx)-repmat(signalX',Nx,1)).^2;
        distance_signalYY = (repmat(signalY,1,Ny)-repmat(signalY',Ny,1)).^2;
        distance_signalXY = (repmat(signalX,1,Ny)-repmat(signalY',Nx,1)).^2;
        
        distance_center_faceXX = zeros(Nx);
        distance_center_faceYY=zeros(Ny);
        distance_center_faceXY=zeros(Nx,Ny);
                
        scp_unit_normalsXX = zeros(Nx);
        scp_unit_normalsYY = zeros(Ny);        
        scp_unit_normalsXY = zeros(Nx,Ny);
    
        for l=1:d
            distance_center_faceXX = distance_center_faceXX+(repmat(center_faceX(:,l),1,Nx)-repmat(center_faceX(:,l)',Nx,1)).^2;
            distance_center_faceYY = distance_center_faceYY+(repmat(center_faceY(:,l),1,Ny)-repmat(center_faceY(:,l)',Ny,1)).^2;
            distance_center_faceXY = distance_center_faceXY+(repmat(center_faceX(:,l),1,Ny)-repmat(center_faceY(:,l)',Nx,1)).^2;

            scp_unit_normalsXX = scp_unit_normalsXX + (repmat(unit_normalsX(:,l),1,Nx).*repmat(unit_normalsX(:,l)',Nx,1));
            scp_unit_normalsYY = scp_unit_normalsYY + (repmat(unit_normalsY(:,l),1,Ny).*repmat(unit_normalsY(:,l)',Ny,1));
            scp_unit_normalsXY = scp_unit_normalsXY + (repmat(unit_normalsX(:,l),1,Ny).*repmat(unit_normalsY(:,l)',Nx,1));
        end
        
        % Geometric kernel      
        Kernel_geomXX = radial_function_geom(distance_center_faceXX,0,objfun);
        Kernel_geomYY = radial_function_geom(distance_center_faceYY,0,objfun);
        Kernel_geomXY = radial_function_geom(distance_center_faceXY,0,objfun);
        
        % Signal kernel
        Kernel_signalXX = radial_function_signal(distance_signalXX,0,objfun);
        Kernel_signalYY = radial_function_signal(distance_signalYY,0,objfun);
        Kernel_signalXY = radial_function_signal(distance_signalXY,0,objfun);

        % tangent space kernel
        Kernel_tanXX = radial_function_sphere(scp_unit_normalsXX,0,objfun);
        Kernel_tanYY =  radial_function_sphere(scp_unit_normalsYY,0,objfun);
        Kernel_tanXY = radial_function_sphere(scp_unit_normalsXY,0,objfun);
        
        % Area 
        totalWX=norm_normalsX.*WX;
        totalWY=norm_normalsY.*WY;   
        AreaXX = (totalWX * totalWX');
        AreaYY = (totalWY * totalWY');
        AreaXY = (totalWX * totalWY');
        
        % norm(x)=
        PXX = sum(sum(AreaXX .* Kernel_geomXX .* Kernel_signalXX .* Kernel_tanXX ));
        
        % morm(y) =
        PYY =sum(sum(AreaYY .* Kernel_geomYY .* Kernel_signalYY .* Kernel_tanYY));
        
        %prs(x,y) =
        PXY = sum(sum(AreaXY .* Kernel_geomXY .* Kernel_signalXY .* Kernel_tanXY));
        
%end         

 g= PXX + PYY - 2* PXY;
end

%------------------------------------%
% equivalent code with scalarProduct %
%------------------------------------%


%function g = shape_Kernel_distance(fs1,fs2,objfun)
% % 
% %Possible method are 'cuda' or 'matlab'.

% % discretize the fshapes
%[center_faceX,signalX,normalsX]=fcatoms(fs1.x,fs1.f,fs1.G,objfun.signal_type);
%[center_faceY,signalY,normalsY]=fcatoms(fs2.x,fs2.f,fs2.G,objfun.data_signal_type);

% % compute the morm squared
%PXX = fshape_kernel_scalarProduct(center_faceX,signalX,normalsX,center_faceX,signalX,normalsX,objfun);
%PYY = fshape_kernel_scalarProduct(center_faceY,signalY,normalsY,center_faceY,signalY,normalsY,objfun);
%PXY = fshape_kernel_scalarProduct(center_faceX,signalX,normalsX,center_faceY,signalY,normalsY,objfun);

%g= PXX + PYY - 2* PXY;

%end
