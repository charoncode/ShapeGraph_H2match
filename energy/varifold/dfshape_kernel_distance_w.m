function [dxg,dfg,drhog]= dfshape_kernel_distance_w(fs1,fs2,objfun)
% DFSHAPE_KERNEL_DISTANCE(templatefinal,target,objfun) computes the derivatives of the kernel based
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
%   dxg : a matrix of the size of templatefinal.x.
%   dfg : a matrix of the size of templatefinal.f.
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2016)


%------------------------%
% Start the chain's rule %
%------------------------%

d=size(fs1.x,2);
m=size(fs1.G,2)-1;

% discretize the fshapes
[center_faceX,signalX,normalsX]=fcatoms(fs1.x,fs1.f,fs1.G,objfun.signal_type);
[center_faceY,signalY,normalsY]=fcatoms(fs2.x,fs2.f,fs2.G,objfun.data_signal_type);

WX=1+fs1.rho;
WY=1+fs2.rho;

% Filter out degenerate faces
% M_areaX=sqrt(sum(normalsX.^2,2));
% ind=find(M_areaX>1e-3*mean(M_areaX));
% center_faceX=center_faceX(ind,:);
% signalX=signalX(ind,:);
% normalsX=normalsX(ind,:);

if min(m,d-m) ==0 % for points clouds or simplexes dim or codim == 0 : some simplifications occurs
    
    [Dcenter_faceXg,DnormalsXg,DsignalXg,drhog] = dptcloud_or_simplexes_kernel_distance(center_faceX,signalX,normalsX,WX,center_faceY,signalY,normalsY,WY,objfun.kernel_distance);
    
elseif min(m,d-m) ==1 % for curves or surface dim or codim ==1;
    
    % add a correcting weight if the target has missing cells
    if isfield(objfun,'weight_missing_data')
        normalsY = normalsY / objfun.weight_missing_data;
    end
    [Dcenter_faceXg,DnormalsXg,DsignalXg,drhog] = dcurves_or_surfaces_kernel_distance(center_faceX,signalX,normalsX,WX,center_faceY,signalY,normalsY,WY,objfun.kernel_distance);
    
end

%-------------------------%
% End of the chain's rule %
%-------------------------%

[dxg,dfg]= chain_rule(fs1.x,fs1.G,Dcenter_faceXg,DnormalsXg,DsignalXg,objfun.signal_type);

end






function [Dcenter_faceXg,DnormalsXg,DsignalXg,drhog]=dptcloud_or_simplexes_kernel_distance(center_faceX,signalX,normalsX,WX,center_faceY,signalY,normalsY,WY,objfun)

% switch objfun.method
%         case 'cuda'  % use cuda to speedup the computation
%     
%      %--------------------------------%
%      % Compute derivative wrt signals %
%      %--------------------------------%
%     	 eval(['dffshape_scp=@dffsimplex_scp_Gpu_',lower(objfun.kernel_geom),lower(objfun.kernel_signal),';']);
%     
%             DsignalXg = 2 *(dffshape_scp(center_faceX',center_faceX',signalX',signalX',normalsX',normalsX',objfun.kernel_size_geom,objfun.kernel_size_signal)...
%                 - dffshape_scp(center_faceX',center_faceY',signalX',signalY',normalsX',normalsY',objfun.kernel_size_geom,objfun.kernel_size_signal))';
%     
%     
%      %--------------------------------%
%      % Compute derivative wrt normals %
%      %--------------------------------%
%     	 eval(['dXifshape_scp=@dXifsimplex_scp_Gpu_',lower(objfun.kernel_geom),lower(objfun.kernel_signal),';']);
%             DnormalsXg = 2 * (dXifshape_scp(center_faceX',center_faceX',signalX',signalX',normalsX',normalsX',objfun.kernel_size_geom,objfun.kernel_size_signal)...
%                 - dXifshape_scp(center_faceX',center_faceY',signalX',signalY',normalsX',normalsY',objfun.kernel_size_geom,objfun.kernel_size_signal))';
%     
%      %------------------------------------%
%      % Compute derivative wrt center_face %
%      %------------------------------------%
%     	 eval(['dXfshape_scp=@dXfsimplex_scp_Gpu_',lower(objfun.kernel_geom),lower(objfun.kernel_signal),';']);
%             Dcenter_faceXg = 2 * (dXfshape_scp(center_faceX',center_faceX',signalX',signalX',normalsX',normalsX',objfun.kernel_size_geom,objfun.kernel_size_signal)...
%                 - dXfshape_scp(center_faceX',center_faceY',signalX',signalY',normalsX',normalsY',objfun.kernel_size_geom,objfun.kernel_size_signal))';
%    
%     otherwise
        
        [Tx,d]=size(center_faceX);
        Ty=size(center_faceY,1);
        
        % compute squared distances and angles
        diff_signalXX = (repmat(signalX,1,Tx)-repmat(signalX',Tx,1));
        diff_signalXY = (repmat(signalX,1,Ty)-repmat(signalY',Tx,1));
        
        distance2_center_faceXX = zeros(Tx);
        distance2_center_faceXY=zeros(Tx,Ty);
        
        for l=1:d
            distance2_center_faceXX = distance2_center_faceXX+(repmat(center_faceX(:,l),1,Tx)-repmat(center_faceX(:,l)',Tx,1)).^2;
            distance2_center_faceXY = distance2_center_faceXY+(repmat(center_faceX(:,l),1,Ty)-repmat(center_faceY(:,l)',Tx,1)).^2;
        end
        
        % compute  Geometric kernel
        Kernel_geomXX  = radial_function_geom(distance2_center_faceXX,0,objfun);
        Kernel_geomXY  = radial_function_geom(distance2_center_faceXY,0,objfun);
        dKernel_geomXX = radial_function_geom(distance2_center_faceXX,1,objfun);
        dKernel_geomXY = radial_function_geom(distance2_center_faceXY,1,objfun);
        
        % compute Signal kernel
        Kernel_signalXX  = radial_function_signal(diff_signalXX.^2,0,objfun);
        Kernel_signalXY  = radial_function_signal(diff_signalXY.^2,0,objfun);
        dKernel_signalXX = radial_function_signal(diff_signalXX.^2,1,objfun);
        dKernel_signalXY = radial_function_signal(diff_signalXY.^2,1,objfun);
        
        % compute Area
        AreaXX = (normalsX * normalsX');
        AreaXY = (normalsX * normalsY');
        totalWX=normalsX.*repmat(WX,1,size(normalsX,2));
        totalWY=normalsY.*repmat(WY,1,size(normalsY,2));
        WXX = (totalWX * totalWX');
        WXY = (totalWX * totalWY');
        
        %------------------------------------%
        % Compute derivative wrt center_face %
        %------------------------------------%
        DXX  = WXX .* dKernel_geomXX .* Kernel_signalXX;
        DXY  = WXY .* dKernel_geomXY .* Kernel_signalXY;
        
        Dcenter_faceXg=zeros(Tx,d);
        for l=1:d
            Dcenter_faceXg(:,l)=4*(sum(( DXX .*(repmat(center_faceX(:,l),1,Tx)-repmat(center_faceX(:,l)',Tx,1))),2)...
                -sum(DXY .* (repmat(center_faceX(:,l),1,Ty)-repmat(center_faceY(:,l)',Tx,1)),2)); % scalar kernel case
        end
        
        %--------------------------------%
        % Compute derivative wrt normals %
        %--------------------------------%
        mXX = Kernel_geomXX .* Kernel_signalXX.*(WX*WX');
        mXY = Kernel_geomXY .* Kernel_signalXY.*(WX*WY');
        
        DnormalsXg = 2*(  mXX * normalsX - mXY * normalsY );
        
        %--------------------------------%
        % Compute derivative wrt signals %
        %--------------------------------%
        DXX  = WXX .* Kernel_geomXX .* dKernel_signalXX .* diff_signalXX;
        DXY  = WXY .* Kernel_geomXY .* dKernel_signalXY .* diff_signalXY;
        
        DsignalXg=4*(sum(DXX,2)-sum(DXY,2));
        
        %--------------------------------%
        % Compute derivative wrt var weights %
        %--------------------------------%        
        
        drhog=2*(sum(AreaXX .* Kernel_geomXX .* Kernel_signalXX .*repmat(WX,1,Tx),2) - sum(AreaXY .* Kernel_geomXY .* Kernel_signalXY .*repmat(WY',Ty,1),2));
%end

end

function [Dcenter_faceXg,DnormalsXg,DsignalXg,drhog]=dcurves_or_surfaces_kernel_distance(center_faceX,signalX,normalsX,WX,center_faceY,signalY,normalsY,WY,objfun)

% switch objfun.method
%     case 'cuda'  % use cuda to speedup the computation
%         
%         %--------------------------------%
%         % Compute derivative wrt signals %
%         %--------------------------------%
%         eval(['dffshape_scp=@dffshape_scp_Gpu_',lower(objfun.kernel_geom),lower(objfun.kernel_signal),lower(objfun.kernel_grass),';']);
%         
%         DsignalXg = 2 *(dffshape_scp(center_faceX',center_faceX',signalX',signalX',normalsX',normalsX',objfun.kernel_size_geom,objfun.kernel_size_signal,objfun.kernel_size_grass)...
%             - dffshape_scp(center_faceX',center_faceY',signalX',signalY',normalsX',normalsY',objfun.kernel_size_geom,objfun.kernel_size_signal,objfun.kernel_size_grass))';
%         
%         
%         %--------------------------------%
%         % Compute derivative wrt normals %
%         %--------------------------------%
%         eval(['dXifshape_scp=@dXifshape_scp_Gpu_',lower(objfun.kernel_geom),lower(objfun.kernel_signal),lower(objfun.kernel_grass),';']);
%         DnormalsXg = 2 * (dXifshape_scp(center_faceX',center_faceX',signalX',signalX',normalsX',normalsX',objfun.kernel_size_geom,objfun.kernel_size_signal,objfun.kernel_size_grass)...
%             - dXifshape_scp(center_faceX',center_faceY',signalX',signalY',normalsX',normalsY',objfun.kernel_size_geom,objfun.kernel_size_signal,objfun.kernel_size_grass))';
%         
%         %------------------------------------%
%         % Compute derivative wrt center_face %
%         %------------------------------------%
%         eval(['dXfshape_scp=@dXfshape_scp_Gpu_',lower(objfun.kernel_geom),lower(objfun.kernel_signal),lower(objfun.kernel_grass),';']);
%         Dcenter_faceXg = 2 * (dXfshape_scp(center_faceX',center_faceX',signalX',signalX',normalsX',normalsX',objfun.kernel_size_geom,objfun.kernel_size_signal,objfun.kernel_size_grass)...
%             - dXfshape_scp(center_faceX',center_faceY',signalX',signalY',normalsX',normalsY',objfun.kernel_size_geom,objfun.kernel_size_signal,objfun.kernel_size_grass))';
%         
%     otherwise
        
        [Tx,d]=size(center_faceX);
        Ty=size(center_faceY,1);
        
        % Compute unit normals
        norm_normalsX = sqrt(sum(normalsX .^2,2));
        norm_normalsY = sqrt(sum(normalsY .^2,2));
        
        unit_normalsX = normalsX ./  repmat(norm_normalsX,1,size(normalsX,2));
        unit_normalsY = normalsY ./  repmat(norm_normalsY,1,size(normalsY,2));
        
        % compute squared distances and angles
        diff_signalXX = (repmat(signalX,1,Tx)-repmat(signalX',Tx,1));
        diff_signalXY = (repmat(signalX,1,Ty)-repmat(signalY',Tx,1));
        
        distance2_center_faceXX = zeros(Tx);
        distance2_center_faceXY=zeros(Tx,Ty);
        
        scp_unit_normalsXX = zeros(Tx);
        scp_unit_normalsXY = zeros(Tx,Ty);
        
        for l=1:d
            distance2_center_faceXX = distance2_center_faceXX+(repmat(center_faceX(:,l),1,Tx)-repmat(center_faceX(:,l)',Tx,1)).^2;
            distance2_center_faceXY = distance2_center_faceXY+(repmat(center_faceX(:,l),1,Ty)-repmat(center_faceY(:,l)',Tx,1)).^2;
            
            scp_unit_normalsXX = scp_unit_normalsXX + (repmat(unit_normalsX(:,l),1,Tx).*repmat(unit_normalsX(:,l)',Tx,1));
            scp_unit_normalsXY = scp_unit_normalsXY + (repmat(unit_normalsX(:,l),1,Ty).*repmat(unit_normalsY(:,l)',Tx,1));
        end
        
        % compute  Geometric kernel
        Kernel_geomXX  = radial_function_geom(distance2_center_faceXX,0,objfun);
        Kernel_geomXY  = radial_function_geom(distance2_center_faceXY,0,objfun);
        dKernel_geomXX = radial_function_geom(distance2_center_faceXX,1,objfun);
        dKernel_geomXY = radial_function_geom(distance2_center_faceXY,1,objfun);
        
        % compute Signal kernel
        Kernel_signalXX  = radial_function_signal(diff_signalXX.^2,0,objfun);
        Kernel_signalXY  = radial_function_signal(diff_signalXY.^2,0,objfun);
        dKernel_signalXX = radial_function_signal(diff_signalXX.^2,1,objfun);
        dKernel_signalXY = radial_function_signal(diff_signalXY.^2,1,objfun);
        
        % compute tangent space kernel
        Kernel_tanXX  = radial_function_sphere(scp_unit_normalsXX,0,objfun);
        Kernel_tanXY  = radial_function_sphere(scp_unit_normalsXY,0,objfun);
        dKernel_tanXX = radial_function_sphere(scp_unit_normalsXX,1,objfun);
        dKernel_tanXY = radial_function_sphere(scp_unit_normalsXY,1,objfun);
        
        % compute Area
        AreaXX = (normalsX * normalsX');
        AreaXY = (normalsX * normalsY');
        
        totalWX=norm_normalsX.*WX;
        totalWY=norm_normalsY.*WY;   
        WXX = (totalWX * totalWX');
        WXY = (totalWX * totalWY');
        
        %------------------------------------%
        % Compute derivative wrt center_face %
        %------------------------------------%
        DXX  = WXX .* dKernel_geomXX .* Kernel_signalXX;
        DXY  = WXY .* dKernel_geomXY .* Kernel_signalXY;
        
        Dcenter_faceXg=zeros(Tx,d);
        for l=1:d
            Dcenter_faceXg(:,l)=4*(sum(( DXX .*(repmat(center_faceX(:,l),1,Tx)-repmat(center_faceX(:,l)',Tx,1))),2)...
                -sum(DXY .* (repmat(center_faceX(:,l),1,Ty)-repmat(center_faceY(:,l)',Tx,1)),2)); % scalar kernel case
        end
        
        %--------------------------------%
        % Compute derivative wrt normals %
        %--------------------------------%
        MXX = Kernel_geomXX .* Kernel_signalXX .* dKernel_tanXX.*(WX*WX');
        MXY = Kernel_geomXY .* Kernel_signalXY .* dKernel_tanXY.*(WX*WY');
        mXX = Kernel_geomXX .* Kernel_signalXX .* Kernel_tanXX.*(WX*WX');
        mXY = Kernel_geomXY .* Kernel_signalXY .* Kernel_tanXY.*(WX*WY');
        
        DnormalsXg = 2*(  repmat(mXX * norm_normalsX  - (MXX .* scp_unit_normalsXX  ) * norm_normalsX,1,d) .* unit_normalsX ...
            + MXX * normalsX ...
            - repmat(mXY * norm_normalsY  - (MXY .* scp_unit_normalsXY  ) * norm_normalsY,1,d) .* unit_normalsX...
            - MXY * normalsY);
        
        %--------------------------------%
        % Compute derivative wrt signals %
        %--------------------------------%
        DXX  = WXX .* Kernel_geomXX .* Kernel_tanXX .* dKernel_signalXX .* diff_signalXX;
        DXY  = WXY .* Kernel_geomXY .* Kernel_tanXY .* dKernel_signalXY .* diff_signalXY;
        
        DsignalXg=4*(sum(DXX,2)-sum(DXY,2));
        
        drhog=2*(sum(AreaXX .* Kernel_geomXX .* Kernel_tanXX.* Kernel_signalXX .*repmat(WX,1,Tx),2) - sum(AreaXY .* Kernel_geomXY.* Kernel_tanXY .* Kernel_signalXY .*repmat(WY',Tx,1),2));
%end

end









