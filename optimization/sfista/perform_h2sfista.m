% perform_h2sfista.m: function that implements the SFISTA algorithm to find the optimal geodesic path and transformed 
% source and/or target edge weights to minimize the matching functional for performing shape graph registration using 
% second order elastic Sobolev metrics (H2 metrics.
%
% Input:
%   initvec   : initializations for the SFISTA procedure
%   optim     : structure containing optimization information
%
% Output:
%   X         : long vector containing coordinates of the path of shape graphs (curvePath), the source and/or
%               target edge weights (rho0 and/or rho1), rotation angle (beta), translation vector (v)
%               and scaling factor (kappa). [(d*N*(Nt-1)+|E0|+|E1|+d+2) x 1 array]
%   summary   : structure containing costs and information from the optimization process
%
%   Note: The first d*N*(Nt-1) entries of X correspond to the coordinates of the path (of component curves) of the 
%   transformed source shape graph, excluding the source (curvePath = [c^k(t)])), the next |E0| and/or |E1| entries
%   correspond to the transformed source edge weights rho0 and/or target edge weights rho1 respectively, and the last 
%   d+2 entries correspond to the rotation angle (beta), translation vector (v) and scaling factor (kappa), in that order.
%   That is, we have: 
%   X = [c^1(2)(:,1);...; c^K(Nt)(:,1); c^1(2)(:,2);...; c^K(Nt)(:,2); rho0; rho1; beta; v(:); kappa]

function [X, summary] = perform_h2sfista(initvec,objfun,optim)

global template data gamma

% Rename stopping tolerance
tol = optim.sfista.tol;

% Dimensions
n = length(initvec); 
[E0,~] = size(template.G); 
[E1,~] = size(data.G); 
[~, d] = size(template.x);
numControlPts = 0;  % number of spatial control points for transformed source
for k = 1:length(template.connComp)
    numControlPts = numControlPts + size(template.connComp{k}.cPts, 1);    
end

% Spline parameters
splineData = objfun.splineData;
Nt = splineData{1}.Nt;

% Extract rotation angle, translation vector, and scaling factor and specify SFISTA initializations
beta_tp1 = initvec(n-d-1);  % beta, rotation angle
v_tp1 = initvec(n-d:n-1);  % v, translation vector
kappa_tp1 = initvec(n);  % kappa, scaling factor

% Extract path of shape graphs (curvePath) and specify SFISTA initializations
iter = 0; 
curvePath_tp1 = initvec(1:d*numControlPts*(Nt-1)); 
curvePath_t = curvePath_tp1 + 2*tol;

% Extract transformed source and/or target edge weights and specify SFISTA initializations
switch lower(objfun.edge_weight)
    
    case 'both'  % if optimizing over both transformed source and target edge weights  
        gamma = [optim.sfista.gamma_source, optim.sfista.gamma_target];  % define as global variable
        rho0_tp1 = initvec(d*numControlPts*(Nt-1)+1:d*numControlPts*(Nt-1)+E0); 
        rho0_t = rho0_tp1 + 2*tol;
        rho1_tp1 = initvec(d*numControlPts*(Nt-1)+E0+1:d*numControlPts*(Nt-1)+E0+E1);
        rho1_t = rho1_tp1 + 2*tol;
    
    case 'source'  % if optimizing over transformed source edge weights only
        gamma = optim.sfista.gamma_source;
        rho0_tp1 = initvec(d*numControlPts*(Nt-1)+1:d*numControlPts*(Nt-1)+E0); 
        rho0_t = rho0_tp1 + 2*tol;
        
    case 'target'  % if optimizing over target edge weights only
        gamma = optim.sfista.gamma_target;
        rho1_tp1 = initvec(d*numControlPts*(Nt-1)+1:d*numControlPts*(Nt-1)+E1); 
        rho1_t = rho1_tp1 + 2*tol;
        
end

% Run SFISTA to optimize over transformed source and/or target edge weights
switch lower(objfun.edge_weight)

    case 'both'  % if optimizing over both transformed source and target edge weights  
        
        % Until convergence criteria is met
        while (norm(curvePath_tp1 - curvePath_t) > tol) && (norm(rho0_tp1 - rho0_t) > tol) &&...
          (norm(rho1_tp1 - rho1_t) > tol) && (iter <= optim.sfista.maxit_gamma)
    
            % Update iteration count
            iter = iter + 1;
    
            % Update optimization variables and reorganize into long vector
            curvePath_t = curvePath_tp1;
            rho0_t = rho0_tp1;
            rho1_t = rho1_tp1;
            beta_t = beta_tp1;
            v_t = v_tp1;
            kappa_t = kappa_tp1;
            x_t = [curvePath_t; rho0_t; rho1_t; beta_t; v_t; kappa_t];
    
            % Increase SFISTA gain parameter according to update rule
            gamma = (1 + sqrt(1 + 4*gamma.^2))./2;
    
            % Call BFGS to minimize the smoothed matching functional wrt curvePath, rho0, rho1, beta, v, kappa
            [x_tp1, summary] = perform_bfgs(optim.sfista.enr_name, x_t, optim.sfista);
    
            % Rearrange x_tp1 for updates
            curvePath_tp1 = x_tp1(1:d*numControlPts*(Nt-1));
            rho0_tp1 = x_tp1(d*numControlPts*(Nt-1)+1:d*numControlPts*(Nt-1)+E0);
            rho1_tp1 = x_tp1(d*numControlPts*(Nt-1)+E0+1:d*numControlPts*(Nt-1)+E0+E1);
            beta_tp1 = x_tp1(end-d-1);
            v_tp1 = x_tp1(end-d:end-1);
            kappa_tp1 = x_tp1(end);
    
        end
      

   case 'source'  % if optimizing over transformed source edge weights only

       % Until convergence criteria is met
       while (norm(curvePath_tp1 - curvePath_t) > tol) && (norm(rho0_tp1 - rho0_t) > tol)  &&...
              (iter <= optim.sfista.maxit_gamma)
    
            % Update iteration count
            iter = iter + 1;
    
            % Update transformed source and edge weights, and reorganize into long vector
            curvePath_t = curvePath_tp1;
            rho0_t = rho0_tp1;
            beta_t = beta_tp1;
            v_t = v_tp1;
            kappa_t = kappa_tp1;
            x_t = [curvePath_t; rho0_t; beta_t; v_t; kappa_t];
    
            % Increase SFISTA gain parameter according to update rule
            gamma = (1 + sqrt(1 + 4*gamma.^2))./2;
    
            % Call BFGS to minimize the smoothed matching functional wrt curvePath, rho0, beta, v, kappa
            [x_tp1, summary] = perform_bfgs(optim.sfista.enr_name, x_t, optim.sfista);
    
            % Rearrange x_tp1 for updates
            curvePath_tp1 = x_tp1(1:d*numControlPts*(Nt-1));
            rho0_tp1 = x_tp1(d*numControlPts*(Nt-1)+1:d*numControlPts*(Nt-1)+E0);
            beta_tp1 = x_tp1(end-d-1);
            v_tp1 = x_tp1(end-d:end-1);
            kappa_tp1 = x_tp1(end);
        
       end
       
        
    case 'target'  % if optimizing over target edge weights only
        
       % Until convergence criteria is met
       while (norm(curvePath_tp1 - curvePath_t) > tol) && (norm(rho1_tp1 - rho1_t) > tol)  &&...
              (iter <= optim.sfista.maxit_gamma)
    
            % Update iteration count
            iter = iter + 1;
    
            % Update transformed source and edge weights, and reorganize into long vector
            curvePath_t = curvePath_tp1;
            rho1_t = rho1_tp1;
            beta_t = beta_tp1;
            v_t = v_tp1;
            kappa_t = kappa_tp1;
            x_t = [curvePath_t; rho1_t; beta_t; v_t; kappa_t];
    
            % Increase SFISTA gain parameter according to update rule
            gamma = (1 + sqrt(1 + 4*gamma.^2))./2;
    
            % Call BFGS to minimize the smoothed matching functional wrt curvePath, rho1, beta, v, kappa
            [x_tp1, summary] = perform_bfgs(optim.sfista.enr_name, x_t, optim.sfista);
    
            % Rearrange x_tp1 for updates
            curvePath_tp1 = x_tp1(1:d*numControlPts*(Nt-1));
            rho1_tp1 = x_tp1(d*numControlPts*(Nt-1)+1:d*numControlPts*(Nt-1)+E1);
            beta_tp1 = x_tp1(end-d-1);
            v_tp1 = x_tp1(end-d:end-1);
            kappa_tp1 = x_tp1(end);
        
       end
       
       
    case 'fixed_weights'
        
       % Until convergence criteria is met
       while (norm(curvePath_tp1 - curvePath_t) > tol) && (iter <= optim.sfista.maxit_gamma)
           
           % Update iteration count
            iter = iter + 1;
    
            % Update transformed source, and reorganize into long vector
            curvePath_t = curvePath_tp1;
            beta_t = beta_tp1;
            v_t = v_tp1;
            kappa_t = kappa_tp1;
            x_t = [curvePath_t; beta_t; v_t; kappa_t];
    
            % Increase SFISTA gain parameter according to update rule
            gamma = (1 + sqrt(1 + 4*gamma.^2))./2;
    
            % Call BFGS to minimize the smoothed matching functional wrt curvePath, rho1, beta, v, kappa
            [x_tp1, summary] = perform_bfgs(optim.sfista.enr_name, x_t, optim.sfista);
    
            % Rearrange x_tp1 for updates
            curvePath_tp1 = x_tp1(1:d*numControlPts*(Nt-1));
            beta_tp1 = x_tp1(end-d-1);
            v_tp1 = x_tp1(end-d:end-1);
            kappa_tp1 = x_tp1(end);
            
       end
    
end 

% Return output
X = x_tp1;

end

