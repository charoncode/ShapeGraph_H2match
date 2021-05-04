% enr_shapegraph_h2match_sfista.m: function that evaluates:
% a) the discretized matching functional for shape graph registration using second order elastic Sobolev
%    metrics (H2 metrics) - to be used when optimizing edge weights using SFISTA.     
% b) the differential of the matching functional wrt the coordinates of the path (of the component curves) of the 
%    transformed source shape graph (curvePath = [c^k(t)] for k = 1,...,K and t = 2,...,Nt), wrt the transformed source 
%    and/or transformed target edge weights (rho0 and/or rho1), and wrt to rotations (beta), translations (v) and scalings (kappa).                  
%
% Input:
%   x   : long vector containing coordinates of the path of curves (of each component) of the transformed source (curvePath),
%         the transformed source and/or target edge weights (rho0 and/or rho1), rotation angle (beta), translation vector (v)
%         and scaling factor (kappa). [(d*N*(Nt-1)+|E0|+|E1|+d+2) x 1 array]
%   Note: The first d*N*(Nt-1) entries of x correspond to the coordinates of the path of curves 
%         excluding the source shape graph (curvePath = [c^k(t)])), the next |E0| and/or |E1| entries
%         correspond to rho0 and/or rho1 respectively, and the last d+2 entries correspond to the 
%         rotation angle (beta), translation vector (v) and scaling factor (kappa), in that order.
%         That is, we have: 
%         x = [c^1(2)(:,1);...; c^K(Nt)(:,1); c^1(2)(:,2);...; c^K(Nt)(:,2); rho0; rho1; beta; v; kappa]
% 
% Output:
%   ENR : value of the matching functional evaluated at current iterates of curvePath,
%         rho0 and/or rho1, beta, v and kappa. [scalar]
%   dENR: differential of the energy functional evaluated at curvePath, rho0 and/or rho1,
%         beta, v and kappa [(d*N*(Nt-1)+|E0|+|E1|+d+2)x1 array]
%   Note: The first d*N*(Nt-1) entries of dENR correspond to the differential wrt to curvePath,
%         the next |E0| and/or |E1| entries correspond to gradients wrt to rho0 and/or rho1 
%         respectively, and the last d+2 entries correspond to differentials wrt to beta, v
%         and kappa, in that order.
%
% Note: This is a simple wrapper function that can be called by Hanso bfgs.


function [ENR, dENR] = enr_shapegraph_h2match_sfista(x)


% Note: template = source (c0), and data = target (c1)
global template data objfunc optimc gamma gradavg


%-----------------------%
%       PARAMETERS      %
%-----------------------%

% Dimensions
n = length(x); 
[E0, ~] = size(template.G); 
[E1, ~] = size(data.G);
[nTarget, d] = size(data.x);  % number of vertices for target
numControlPts = 0;  % number of spatial control points for transformed source
for k = 1:length(template.connComp)
    numControlPts = numControlPts + size(template.connComp{k}.cPts, 1);    
end

% Extract structures containing the matching functional and optimization parameters
objfun = objfunc; 
optim = optimc; 

% Spline parameters
splineData = objfun.splineData;
Nt = splineData{1}.Nt;

% Varifold parameters
lambdaVar = objfun.kernel_distance.lambda_var;

% TV norm parameters
lambdaTvSource = objfun.tv.lambda_tv_source;
lambdaTvTarget = objfun.tv.lambda_tv_target;

% Double-well penalty parameters
lambdaPenSource = objfun.penalty.lambda_pen_source;
lambdaPenTarget = objfun.penalty.lambda_pen_target;
m0 = objfun.penalty.pen_min1;
m1 = objfun.penalty.pen_min2;

% Optimization parameters
options = optim.options;
optRot = options.optRot;
optTra = options.optTra;
optScal = options.optScal;


%--------------------------------------%
%      EXTRACT & UPDATE VARIABLES      %
%--------------------------------------%

% Extract rotation angle, translation vector, and scaling factor
beta = x(n-d-1);  % beta, rotation angle
v = x(n-d:n-1);  % v, translation vector
kappa = x(n);  % kappa, scaling factor

% Extract path of shape graphs across all components
allPaths = reshape(x(1:d*numControlPts*(Nt-1)), [numControlPts*(Nt-1), d]);

% Separate path of shape graphs by component curve
pathComp = cell(length(template.connComp), 1);
nCompPath = 0;

for k = 1:length(template.connComp) 
    
    % Extract spline data for each component curve
    Nk = splineData{k}.N;
    
    % Extract full path of control points for each component curve (including fixed initial curve, i.e, source)
    pathComp{k} = [ template.connComp{k}.cPts ; ...
                    allPaths(nCompPath+1:nCompPath+(Nk*(Nt-1)),:)];
    nCompPath = nCompPath + Nk*(Nt-1);
    
end

% Extract transformed target, and translate, rotate and scale its vertices
transfTarget = data; 
rotation = [cos(beta), sin(beta); -sin(beta), cos(beta)];
if d == 2
    transfTarget.x = kappa * ( (transfTarget.x + (ones(nTarget, 1) * v(:)')) * rotation );
elseif d == 3
    transfTarget.x = kappa * ( (transfTarget.x + (ones(nTarget, 1) * v(:)')) * eye(d) );  % no rotations for 3D curves
end

% Update target edge weights
switch lower(objfun.edge_weight)
    
    case 'both' % if optimizing over both transformed source and target edge weights
        transfTarget.rho = x(d*numControlPts*(Nt-1)+E0+1:d*numControlPts*(Nt-1)+E0+E1);
    
    case 'target' % if optimizing over target edge weights only
        transfTarget.rho = x(d*numControlPts*(Nt-1)+1:d*numControlPts*(Nt-1)+E1);
    
end


%------------------------------------------------------%
%      COMPUTE H2 ENERGY OF PATH AND DIFFERENTIAL      %
%------------------------------------------------------%

% Compute H2 energy of path of transformed source shape graph and its differential
energyPathComp = zeros(length(template.connComp), 1);
dENRh2_curvePathComp = cell(length(template.connComp), 1);
dENRh2_endCurveComp = cell(length(template.connComp), 1);

for k = 1:length(template.connComp)
    
    % Compute H2 energy and differential for each component curve
    [energyPathComp(k), dENRh2_curvePathComp{k}, dENRh2_endCurveComp{k}] = ...
                              h2_energyPath(pathComp{k}, objfun, optim, splineData{k});
                                
end

% Aggregate H2 energies from each component curve
energyPath = sum(energyPathComp);


%-------------------------------------------------------------------------------------%
%      COMPUTE DATA ATTACHMENT TERM, REGULARIZERS, PENALTIES & THEIR DIFFERENTIALS    %
%-------------------------------------------------------------------------------------%

% Define full transformed source by concatenating its component curves
transfSource = template;
transfSource.x = [];
transfSource.rho = [];
transfSource.connComp = cell(length(template.connComp), 1);

B_endCurveCompS = cell(length(template.connComp), 1);  % quadrature weights for each component curve
B_transfSourceS = [];  % block diagonal matrix of quadrature weights for full transformed source
nCompRho = 0;

for k = 1:length(template.connComp) 
    
    % Extract spline parameters
    Nk = splineData{k}.N;
    
    % Extract end curve for each component
    transfSource.connComp{k} = template.connComp{k};
    transfSource.connComp{k}.x = [];
    transfSource.connComp{k}.cPts = pathComp{k}(end-Nk+1:end,:);  % spline control points of end curve
    
    % Extract end curve edge weights
    switch lower(objfun.edge_weight)
        
        % weights placed on edges between (interpolated) vertices of end curve
        case {'both', 'source'}
            
            transfSource.connComp{k}.rho = x(d*numControlPts*(Nt-1)+nCompRho+1:...
                                             d*numControlPts*(Nt-1)+nCompRho+size(template.connComp{k}.G, 1));
            nCompRho = nCompRho + size(template.connComp{k}.G, 1);
    
    end
    
    % Interpolate vertices for each component curve from control points
    B_endCurveCompS{k} = splineData{k}.quadData.B_endCurveS;
    B_transfSourceS = blkdiag(B_transfSourceS, B_endCurveCompS{k});  % interpolation matrix
    transfSource.connComp{k}.x = B_endCurveCompS{k}*transfSource.connComp{k}.cPts;  % vertices on end curve interpolated from spline control points
    
    % Concatenate vertices and edge weights from each component curve
    transfSource.x = [transfSource.x ; transfSource.connComp{k}.x];
    transfSource.rho = [transfSource.rho ; transfSource.connComp{k}.rho];
    
end

% Compute varifold distance between transformed source (c(Nt)) and transformed target (c1)
varDist = matchterm_w(transfSource, transfTarget, objfun);

switch lower(objfun.edge_weight)
    
    case 'both'  % if optimizing over both transformed source and target edge weights
        
        % Rename SFISTA parameter (stored as global variable)
        gamma_source = gamma(1);
        gamma_target = gamma(2);
        
        % Compute differential of varifold distance wrt vertices of transformed source (c(Nt)) and its edge weights (rho0)
        [dVarDist_transfSource, dVarDist_sourceRho] = dmatchterm_w(transfSource, transfTarget, objfun);
        dVarDist_transfSource = gradavg * B_transfSourceS' * dVarDist_transfSource;  % averaging + chain rule (because of spline interpolation) for gradient of varifold distance wrt transformed source
        
        % Compute differential of varifold distance wrt vertices of target (c1) (for rigid motions gradients) and its edge weights (rho1)
        [dVarDist_target, dVarDist_targetRho] = dmatchterm_w(transfTarget, transfSource, objfun);
        
        % Compute Huber function evaluated at transformed source edge weights rho0, and its associated gradient
        v0 = transfSource.D * transfSource.rho;
        [hubSource, dHubSource_SourceRho] = huber(v0, lambdaTvSource, gamma_source);
        
        % Compute Huber function evaluated at target edge weights rho1, and its associated gradient
        v1 = transfTarget.D * transfTarget.rho;
        [hubTarget, dHubTarget_targetRho] = huber(v1, lambdaTvTarget, gamma_target);
        
        % Compute double-well penalty of rho0 and associated differentials
        [penSource, dPenSource_sourceRho, dPenSource_transfSource] = quartic_pen(transfSource, m0, m1, 'polynomial', 'length', 0.5);
        dPenSource_transfSource = gradavg * B_transfSourceS' * dPenSource_transfSource;
    
        % Compute double-well penalty of rho1 and associated differentials
        [penTarget, dPenTarget_targetRho] = quartic_pen(transfTarget, m0, m1, 'polynomial', 'length', 0.5);
    
    case 'source'  % if optimizing over transformed source edge weights only
        
        % Rename SFISTA parameter (stored as global variable)
        gam = gamma;
        
        % Compute differential of varifold distance wrt vertices of transformed source (c(Nt)) and its edge weights (rho0)
        [dVarDist_transfSource, dVarDist_sourceRho] = dmatchterm_w(transfSource, transfTarget, objfun);
        dVarDist_transfSource = gradavg * B_transfSourceS' * dVarDist_transfSource;  % averaging + chain rule (because of spline interpolation) for gradient of varifold distance wrt transformed source
        
        % Compute differential of varifold distance wrt to vertices of target (c1) (for rigid motions gradients)
        [dVarDist_target, ~] = dmatchterm_w(transfTarget, transfSource, objfun);
        
        % Compute Huber function evaluated at transformed source edge weights rho0, and associated gradient
        v0 = transfSource.D * transfSource.rho;
        [hubSource, dHubSource_SourceRho] = huber(v0, lambdaTvSource, gam);
        
        % Compute double-well penalty of rho0 and associated differentials
        [penSource, dPenSource_sourceRho, dPenSource_transfSource] = quartic_pen(transfSource, m0, m1, 'polynomial', 'length', 0.5);
        dPenSource_transfSource = gradavg * B_transfSourceS' * dPenSource_transfSource;
    
    case 'target'  % if optimizing over target edge weights only
        
        % Rename SFISTA parameter (stored as global variable)
        gam = gamma;

        % Compute differential of varifold distance wrt vertices of transformed source (c(Nt))
        [dVarDist_transfSource, ~] = dmatchterm_w(transfSource, transfTarget, objfun);
        dVarDist_transfSource = gradavg * B_transfSourceS' * dVarDist_transfSource; % averaging + chain rule (because of spline interpolation) for gradient of varifold distance wrt transformed source

        % Compute differential of varifold distance wrt vertices of target (c1) (for rigid motions gradients) and its edge weights (rho1)
        [dVarDist_target, dVarDist_targetRho] = dmatchterm_w(transfTarget, transfSource, objfun);
        
        % Compute Huber function evaluated at target edge weights rho1, and associated gradient
        v1 = transfTarget.D * transfTarget.rho;
        [hubTarget, dHubTarget_targetRho] = huber(v1, lambdaTvTarget, gam);

        % Compute double-well penalty of rho1 and associated differentials
        [penTarget, dPenTarget_targetRho] = quartic_pen(transfTarget, m0, m1, 'polynomial', 'length', 0.5);
     
    case 'fixed_weights'
        
        % Compute differential of varifold distance wrt vertices of transformed source (c(Nt))
        [dVarDist_transfSource, ~] = dmatchterm_w(transfSource, transfTarget, objfun);  
        dVarDist_transfSource = gradavg * B_transfSourceS' * dVarDist_transfSource;  % averaging + chain rule (because of splines) for gradient of varifold distance wrt transformed source
        
        % Compute differential of varifold distance wrt vertices of target (c1) (for rigid motions gradients)
        [dVarDist_target, ~] = dmatchterm_w(transfTarget, transfSource, objfun);
        
end


%---------------------------------------------------------------------%
%      COMPUTE GRADIENTS WRT TO ROTATIONS, TRANSLATIONS & SCALINGS    %
%---------------------------------------------------------------------%

% Note: E(rho,beta,v) = varDist(endCurve - kappa * [R_beta * (transf_target + v)])^2 

if optRot && d == 2  % if optimizing over rotations for 2D curves
    rotationDer = [-sin(beta), cos(beta); -cos(beta), -sin(beta)]; 
    dVarDist_beta = kappa * sum(sum(dVarDist_target .* ((transfTarget.x + (ones(nTarget, 1) * v(:)')) * rotationDer)));
else
    dVarDist_beta = 0;  % in particular, no rotations for 3D curves 
end

if optTra % if optimizing over translations
    dVarDist_v = kappa * rotation * sum(dVarDist_target, 1)';
else
    dVarDist_v = zeros(d, 1);
end

if optScal && d == 2  % if optimizing over scalings for 2D curves
    dVarDist_kappa = sum(sum(dVarDist_target .* ((transfTarget.x + (ones(nTarget, 1) * v(:)')) * rotation)));
elseif optScal && d == 3  % if optimizing over scalings for 3D curves
    dVarDist_kappa = sum(sum(dVarDist_target .* ((transfTarget.x + (ones(nTarget, 1) * v(:)')) * eye(d))));  % no rotations for 3D curves
else
    dVarDist_kappa = 0; 
end


%----------------------------------%
%    PERFORM GRADIENT AVERAGING    %
%----------------------------------%

% Note: 
%   - dVarDist_transfSource and dPenSource_transfSource have already been averaged
%   - dENRh2_curvePathComp, dENRh2_endCurve (currently ordered by component) need to be averaged

grad_avg_path = [];  % averaging matrix for path of curves
dENRh2_curvePath_avg = [];  % reordered gradient of H2 energy wrt path of curves
dENRh2_endCurve_avg = [];  % reordered gradient of H2 energy wrt endcurve

% Reorder dENRh2_curvePathComp by time, i.e., as follows: c^1(t),...,c^K(t) for t = 2,3,...,Nt-1
for t = 1:Nt-2   
    for k = 1:length(template.connComp)    
        
        Nk = splineData{k}.N;
        dENRh2_curvePath_avg = [dENRh2_curvePath_avg ; dENRh2_curvePathComp{k}((t-1)*Nk+1:t*Nk,:)]; 
        
    end
    
    grad_avg_path = blkdiag(grad_avg_path, gradavg);  % concatenate gradient averaging matrix for path of curves 
    
end

% Reorder dENRh2_endCurve by time, i.e, as follows: c^1(Nt),...,c^K(Nt)
for k = 1:length(template.connComp)
    
    dENRh2_endCurve_avg = [dENRh2_endCurve_avg ; dENRh2_endCurveComp{k}];
    
end

% Apply averaging matrix to reordered gradients
dENRh2_curvePath_avg = grad_avg_path * dENRh2_curvePath_avg;
dENRh2_endCurve_avg = gradavg * dENRh2_endCurve_avg;

% Re-reorder gradients component by component 
grad_path = cell(length(template.connComp), 1);
grad_endcurve = cell(length(template.connComp), 1);

% Reorder dENRh2_curvePathComp component by component after averaging has been applied
n_endcurve = length(template.cPts);  % total number of control points for endcurve
for t = 1:Nt-2    
    
    ind = 0;  % indexing variable
    
    for k = 1:length(template.connComp)
        
        Nk = splineData{k}.N;
        grad_path{k} = [grad_path{k} ; dENRh2_curvePath_avg((t-1)*n_endcurve+ind+1:(t-1)*n_endcurve+ind+Nk,:)];
        ind = ind + Nk;
        
    end    
end
dENRh2_curvePathComp = grad_path;

% Reorder dENRh2_endCurveComp component by component after averaging has been applied
ind = 0;
for k = 1:length(template.connComp)
    
    Nk = splineData{k}.N;
    grad_endcurve{k} = dENRh2_endCurve_avg(ind+1:ind+Nk,:);
    ind = ind + Nk;
    
end
dENRh2_endCurveComp = grad_endcurve;

    
%-------------------------------------------------------%
%      EVALUATE MATCHING FUNCTIONAL AND ITS GRADIENT    %
%-------------------------------------------------------%

% Toggle for geometric optimization
if optim.options.optGeom 
    u = 1;
else
    u = 0;
end

% Evaluate matching functional and its differential
switch lower(objfun.edge_weight)
    
    case 'both'  % if optimizing over both transformed source and target edge weights
        
        % Evaluate matching functional
        ENR = energyPath + lambdaVar*varDist + ...
                 sum(hubSource) + sum(hubTarget) + ...
                    lambdaPenSource*sum(penSource) + lambdaPenTarget*sum(penTarget);
            
        % Append derivative of matching functional wrt end curves to that of derivative wrt path of component curves
        dENRh2_curvePath = [];
        ind = 0;
        for k = 1:length(template.connComp)    
            Nk = splineData{k}.N;
            dENRh2_curvePath = [dENRh2_curvePath ; dENRh2_curvePathComp{k} ;
                                dENRh2_endCurveComp{k} + lambdaVar*dVarDist_transfSource(ind+1:ind+Nk,:) + ...
                                    lambdaPenSource*dPenSource_transfSource(ind+1:ind+Nk,:)];  % ordered component by component
            ind = ind + Nk;
        end
        
        % Return dENR as a long (d*N*(Nt-1)+E0+E1+d+2)x1 vector with all differentials included
        dENR = [ u*dENRh2_curvePath(:); ...  % wrt path of component curves of shape graph
                 lambdaVar*dVarDist_sourceRho + (transfSource.D)'*dHubSource_SourceRho + ...
                     lambdaPenSource*dPenSource_sourceRho; ...  % wrt to source edge weights
                 lambdaVar*dVarDist_targetRho + (transfTarget.D)'*dHubTarget_targetRho + ...
                     lambdaPenTarget*dPenTarget_targetRho; ...  % wrt to target edge weights
                 lambdaVar*dVarDist_beta; ...  % wrt rotations
                 lambdaVar*dVarDist_v; ...  % wrt translations
                 lambdaVar*dVarDist_kappa ];  % wrt scalings
    
    case 'source'  % if optimizing over transformed source edge weights only
        
        % Evaluate matching functional
        ENR = energyPath + lambdaVar*varDist + ...
                     sum(hubSource) + lambdaPenSource*sum(penSource);
            
        % Append derivative of matching functional wrt end curves to that of derivative wrt path of component curves
        dENRh2_curvePath = [];
        ind = 0;
        for k = 1:length(template.connComp)    
            Nk = splineData{k}.N;
            dENRh2_curvePath = [dENRh2_curvePath ; dENRh2_curvePathComp{k} ;
                                dENRh2_endCurveComp{k} + lambdaVar*dVarDist_transfSource(ind+1:ind+Nk,:) + ...
                                    lambdaPenSource*dPenSource_transfSource(ind+1:ind+Nk,:)];  % ordered component by component
            ind = ind + Nk;
        end
        
        % Return dENR as a long (d*N*(Nt-1)+E0+d+2)x1 vector with all differentials included
        dENR = [ u*dENRh2_curvePath(:); ...  % wrt path of component curves of shape graph
                 lambdaVar*dVarDist_sourceRho + (transfSource.D)'*dHubSource_SourceRho + ...
                     lambdaPenSource*dPenSource_sourceRho; ...  % wrt to source edge weights
                 lambdaVar*dVarDist_beta; ...  % wrt rotations
                 lambdaVar*dVarDist_v; ...  % wrt translations
                 lambdaVar*dVarDist_kappa ];  % wrt scalings
             
    case 'target'  % if optimizing over target edge weights only

        % Evaluate matching functional
        ENR = energyPath + lambdaVar*varDist + ... 
                    sum(hubTarget) + lambdaPenTarget*sum(penTarget);
        
        % Append derivative of matching functional wrt end curves to that of derivative wrt path of component curves
        dENRh2_curvePath = [];
        ind = 0;
        for k = 1:length(template.connComp)    
            Nk = splineData{k}.N;
            dENRh2_curvePath = [dENRh2_curvePath ; dENRh2_curvePathComp{k} ;
                                dENRh2_endCurveComp{k} + lambdaVar*dVarDist_transfSource(ind+1:ind+Nk,:)];  % ordered component by component
            ind = ind + Nk;
        end
        
        % Return dENR as a long (d*N*(Nt-1)+E1+d+2)x1 vector with all differentials included
        dENR = [ u*dENRh2_curvePath(:); ...  % wrt path of component curves of shape graph
                 lambdaVar*dVarDist_targetRho + (transfTarget.D)'*dHubTarget_targetRho + ...
                     lambdaPenTarget*dPenTarget_targetRho; ...  % wrt to target edge weights
                 lambdaVar*dVarDist_beta; ...  % wrt rotations
                 lambdaVar*dVarDist_v; ...  % wrt translations
                 lambdaVar*dVarDist_kappa ];  % wrt scalings
             
    case 'fixed_weights' % fixed weights on transformed source and target, no optimization over weights
        
        % Evaluate matching functional
        ENR = energyPath + lambdaVar*varDist;
        
        % Append derivative of matching functional wrt end curves to that of derivative wrt path of component curves
        dENRh2_curvePath = [];
        ind = 0;
        for k = 1:length(template.connComp)    
            Nk = splineData{k}.N;
            dENRh2_curvePath = [dENRh2_curvePath ; dENRh2_curvePathComp{k} ;
                                dENRh2_endCurveComp{k} + lambdaVar*dVarDist_transfSource(ind+1:ind+Nk,:)];  % ordered component by component
            ind = ind + Nk;
        end
        
        % Return dENR as a long (d*N*(Nt-1)+d+2)x1 vector with all differentials included
        dENR = [ u*dENRh2_curvePath(:); ...  % wrt path of component curves of shape graph
                 lambdaVar*dVarDist_beta; ...  % wrt rotations
                 lambdaVar*dVarDist_v; ...  % wrt translations
                 lambdaVar*dVarDist_kappa ];  % wrt scalings
        
end


end

