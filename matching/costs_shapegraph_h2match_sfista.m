% costs_shapegraph_h2match_sfista.m: function that evaluates the costs in the discretized matching functional for
% shape graph registration using second order elastic Sobolev metrics (H2 metrics) - when SFISTA is used.
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
%   ENR: structure containing values of costs and energies appearing in the discretized matching functional, namely:
%   
%       a) energy_path  : value of the energy of the geodesic path of the transformed source
%       b) varifold_cost: value of the varifold distance between the transformed source c(1) and target c1
%       c) huber_source : huber penalty of transformed source edge weights rho0
%       d) huber_target : huber penalty of target edge weights rho1
%       e) pen_source   : double-well penalty of transformed source edge weights rho0
%       f) pen_target   : double-well penalty of target edge weights rho1
%       g) energy_min   : sum of all components of the matching functional weighted by their coefficients
%       h) opt_path     : optimal geodesic path of (each component curve) of the transformed source (curvePath)
%       i) opt_beta     : optimal rotation angle (beta)
%       j) opt_v        : optimal translation vector (v)
%       h) opt_kappa    : optimal scaling factor (kappa)
%
%   optPath     : structure containing optimal deformation path of transformed source
%   transfSource: structure containing transformed source
%   transfTarget: structure containing transformed target

function [ENR, optPath, transfSource, transfTarget] = costs_shapegraph_h2match_sfista(x)


% Note: template = source (c0), and data = target (c1)
global template data objfunc optimc gamma


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
optPath = cell(length(template.connComp), 1);
nCompPath = 0;

for k = 1:length(template.connComp)
    
    % Extract spline data for each component curve
    Nk = splineData{k}.N;
    
    % Extract full path of control points for each component curve (including fixed initial curve, i.e, source)
    optPath{k} = [ template.connComp{k}.cPts ; ...
                   allPaths(nCompPath+1:nCompPath+Nk*(Nt-1),:)];
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


%--------------------------------------%
%      COMPUTE H2 ENERGY OF PATH       %
%--------------------------------------%

% Compute H2 energy of optimal path of transformed source shape graph
energyPathComp = zeros(length(template.connComp),1);

for k = 1:length(template.connComp)
    
    % Compute H2 energy and differential for each component
    [energyPathComp(k), ~, ~] = h2_energyPath(optPath{k}, objfun, optim, splineData{k});
                                
end

% Aggregate H2 energies from each component
ENR.energy_path = sum(energyPathComp);


%-----------------------------------------------------------------------%
%      COMPUTE ENERGY, DATA ATTACHMENT TERM, REGULARIZERS & PENALTIES   %
%-----------------------------------------------------------------------%

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
    transfSource.connComp{k}.cPts = optPath{k}(end-Nk+1:end,:);  % spline control points of end curve
    
    % Extract end curve edge weights
    switch lower(objfun.edge_weight)
        
        % weights placed on edges between (interpolated) vertices of end curve
        case {'both', 'source'}
            
            transfSource.connComp{k}.rho = x(d*numControlPts*(Nt-1)+nCompRho+1:...
                                             d*numControlPts*(Nt-1)+nCompRho+size(template.connComp{k}.G,1));
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
ENR.varifold_cost = matchterm_w(transfSource, transfTarget, objfun);  

switch lower(objfun.edge_weight)
    
    case 'both'  % if optimizing over both transformed source and target edge weights
        
        % Rename SFISTA gain parameter (stored as global variable)
        gam_source = gamma(1);
        gam_target = gamma(2);
        
        % Compute Huber function of transformed source edge weights rho0
        v0 = transfSource.D * transfSource.rho;
        [hubSource, ~] = huber(v0, lambdaTvSource, gam_source);
        ENR.huber_source = sum(hubSource);
        
        % Compute Huber function of target edge weights rho1
        v1 = transfTarget.D * transfTarget.rho;
        [hubTarget, ~] = huber(v1, lambdaTvTarget, gam_target);
        ENR.huber_target = sum(hubTarget);
        
        % Compute double-well penalty of rho0
        [pen_source, ~, ~] = quartic_pen(transfSource, m0, m1, 'polynomial', 'length', 0.5);
        ENR.pen_source = sum(pen_source);
        
        % Compute double-well penalty of rho1
        [pen_target, ~] = quartic_pen(transfTarget, m0, m1, 'polynomial', 'length', 0.5);
        ENR.pen_target = sum(pen_target);
        
        % Evaluate matching functional
        ENR.energy_min = ENR.energy_path + lambdaVar*ENR.varifold_cost + ...
                 ENR.huber_source + ENR.huber_target + ...
                    lambdaPenSource*ENR.pen_source + lambdaPenTarget*ENR.pen_target;
    
    case 'source'  % if optimizing over transformed source edge weights only
        
        % Rename SFISTA gain parameter (stored as global variable)
        gam_source = gamma;
        
        % Compute Huber function of transformed source edge weights rho0
        v0 = transfSource.D * transfSource.rho;
        [hubSource, ~] = huber(v0, lambdaTvSource, gam_source);
        ENR.huber_source = sum(hubSource);
        
        % Compute double-well penalty of rho0
        [pen_source, ~, ~] = quartic_pen(transfSource, m0, m1, 'polynomial', 'length', 0.5);
        ENR.pen_source = sum(pen_source);
        
        % Evaluate matching functional
        ENR.energy_min = ENR.energy_path + lambdaVar*ENR.varifold_cost + ...
                 ENR.huber_source + lambdaPenSource*ENR.pen_source;
        
    case 'target'  % if optimizing over target edge weights only
        
        % Rename SFISTA gain parameter (stored as global variable)
        gam_target = gamma;
        
        % Compute Huber function of target edge weights rho1
        v1 = transfTarget.D * transfTarget.rho;
        [hubTarget, ~] = huber(v1, lambdaTvTarget, gam_target);
        ENR.huber_target = sum(hubTarget);

        % Compute double-well penalty of rho1
        [pen_target, ~] = quartic_pen(transfTarget, m0, m1, 'polynomial', 'length', 0.5);
        ENR.pen_target = sum(pen_target);
        
        % Evaluate matching functional
        ENR.energy_min = ENR.energy_path + lambdaVar*ENR.varifold_cost + ...
                 ENR.huber_target + lambdaPenTarget*ENR.pen_target;
        
    case 'fixed_weights'  % fixed weights on transformed source and target, no optimization over weights
        
        % Evaluate matching functional
        ENR.energy_min = ENR.energy_path + lambdaVar*ENR.varifold_cost;
     
end


%---------------------------------------------%
%      STORE OPTIMAL RIGID MOTION VARIABLES   %
%---------------------------------------------%

% Store rotation angle, translation vector and scaling factor
ENR.opt_beta = beta;
ENR.opt_v = v;
ENR.opt_kappa = kappa;

end


